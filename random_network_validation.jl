##
using DataFrames
using CSV
using Graphs
using StatsBase
using LinearAlgebra
using CairoMakie
using HypothesisTests
using MultipleTesting
using Pipe
using DelimitedFiles
##
project_name = "TCGA-LUSC"
df_nodes_path = "graph_visualize/" * project_name * "_methylation_tumor_gene_loc_giant.csv"
df_edges_path = "graph_visualize/" * project_name * "_methylation_tumor_gene_giant.csv"

df_nodes = CSV.read(
    df_nodes_path,
    DataFrame,
)

df_edges = CSV.read(
    df_edges_path,
    DataFrame,
)
##
function calculate_C_xy(df_nodes, df_edges; n = 10)
    n_community = length(unique(df_nodes[!, :leiden_membership]))
    connection_matrix = zeros(Int, n_community, n_community)
    n_intercommunity_list = zeros(Int, n_community)
    for row in eachrow(df_edges)
        node1, node2 = row.Target_member, row.Source_member
        connection_matrix[node1, node2] += 1
        connection_matrix[node2, node1] += 1
        if node1 != node2
            n_intercommunity_list[node1] += 1
            n_intercommunity_list[node2] += 1
        end
    end
    for row in eachrow(df_edges)
        node1, node2 = row.Target_member, row.Source_member
        connection_matrix[node1, node2] += 1
        connection_matrix[node2, node1] += 1
        if node1 != node2
            n_intercommunity_list[node1] += 1
            n_intercommunity_list[node2] += 1
        end
    end
    connection_matrix[diagind(connection_matrix)] .= 0
    relative_connection_matrix = connection_matrix ./ n_intercommunity_list
    return (relative_connection_matrix[1:n, 1:n])
end

relative_connection_matrix = calculate_C_xy(df_nodes, df_edges)
##
fig = Figure()
ax = Axis(fig[1, 1])
hm = heatmap!(ax, relative_connection_matrix', colormap = Reverse(:RdYlBu_10))
Colorbar(fig[1, 2], hm)
fig
##

function modified_random_network(df_nodes, df_edges)
    r_network = copy(df_edges[!, [:Target_member, :Source_member]])
    for row in eachrow(r_network)
        if row.Target_member == row.Source_member
            continue
        else
            node1, node2 = sample(1:nrow(df_nodes), 2, replace = false)

            row.Target_member = df_nodes[node1, "leiden_membership"]
            row.Source_member = df_nodes[node2, "leiden_membership"]
        end
    end
    return (r_network)
end
##
n_r_network = 5000

relative_connection_tensor = zeros(10, 10, n_r_network)

for i in 1:n_r_network
    df_edges_r_network = modified_random_network(df_nodes, df_edges)
    relative_connection_tensor[:, :, i] .= calculate_C_xy(df_nodes, df_edges_r_network)
end
##
p_values_matrix = similar(relative_connection_matrix)
for j in 1:size(relative_connection_matrix, 2)
    for i in 1:size(relative_connection_matrix, 1)
        test = OneSampleTTest(relative_connection_tensor[i, j, :], relative_connection_matrix[i, j])
        p_values_matrix[i, j] = pvalue(test)
    end
end

adjusted_p_value = @pipe adjust(PValues(vec(p_values_matrix)), Bonferroni()) |>
                         reshape(_, size(p_values_matrix))

writedlm("cxy_matrix", relative_connection_matrix, ",")
writedlm("cxy_p_value_matrix", adjusted_p_value, ",")
