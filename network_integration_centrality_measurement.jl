##
using Distributed
Distributed.addprocs(20, exeflags = `--project=$(Base.active_project())`)

@everywhere using Graphs
import Graphs.Parallel
@everywhere using MetaGraphs
using GraphIO
using DataFrames
using CSV
using Statistics
using StatsBase
using MetaGraphs
include("leiden_community.jl")
import .JuliaCommunity as juliac
using Pipe
using Random
Random.seed!(1);
##

project_name = "TCGA-LUSC"
tumor_input_network_path = "output_network/output_" * project_name * "_Tumor.tsv"
input_hypermethylation_path = "Result/LUSC/TN_Tumor_vs_Normal/hyper/getPair.hyper.pairs.significant.csv"
input_hypomethylation_path = "Result/LUSC/TN_Tumor_vs_Normal/hypo/getPair.hypo.pairs.significant.csv"
input_gene_loc_path = "output_network/" * project_name * "_gene_location.csv"

tumor_output_gene_loc_path = "graph_visualize/" * project_name * "_methylation_tumor_gene_loc_giant.csv"
tumor_output_network_path = "graph_visualize/" * project_name * "_methylation_tumor_gene_giant.csv"
output_graph = "graph_visualize/" * project_name * "_meth_gene_graph_giant"
##
# Helper function to integrate and perform community identification
function get_df_network(path; threshold = 0.99, delim = "\t")
    df = CSV.read(
        path,
        DataFrame,
        header = ["Source", "Target", "Importance"],
        delim = delim,
    )
    nrow = size(df, 1)
    df = df[1:2:nrow, :]
    limit = quantile(df.Importance, threshold)
    df = df[df.Importance.>limit, :]
    return df
end

function get_df_methylation(; hypermethylation_path, hypomethylation_path)
    df_hypomethyl = CSV.read(
        hypomethylation_path,
        DataFrame,
        header = 1,
        select = ["Probe", "GeneID", "Distance"],
    )
    df_hypomethyl[!, :link_type] .= "hypomethylation"
    df_hypermethyl = CSV.read(
        hypermethylation_path,
        DataFrame,
        header = 1,
        select = ["Probe", "GeneID", "Distance"],
    )
    df_hypermethyl[!, :link_type] .= "hypermethylation"
    append!(df_hypomethyl, df_hypermethyl)
    rename!(df_hypomethyl, Dict(:Probe => :Source, :GeneID => :Target))
    return df_hypomethyl
end

function get_gene_location_dictionary(node_set, df_gene_loc)
    gene_loc_dict = Dict()
    for row in eachrow(df_gene_loc)
        if row.ensembl_gene_id in node_set
            gene_loc_dict[row.ensembl_gene_id] = (row.external_gene_name, row.seqnames, row.start, row.end)
        end
    end
    return gene_loc_dict
end

function add_methylation_node!(gene_location_dict, df_methylation)
    for row in eachrow(df_methylation)
        gene_location_dict[row.Source] = ("methylation", gene_location_dict[row.Target][2], -1, -1)
    end
end

function create_node_set_dataframe(gene_location_dict)
    df = DataFrame(id = String[], label = String[], loc = String[], type = [])
    for (key, value) in gene_location_dict
        if value[1] == "methylation"
            push!(df, (key, key, value[2], "probe"))
        else
            push!(df, (key, value[1], value[2], "gene"))
        end
    end
    return df
end

function create_gene_methylation_network_dataframe(df_gene_network, df_methylation, location_dict)
    #Create gene network dataframe
    df_gene = select(df_gene_network, [:Source, :Target])
    distance = zeros(Union{Missing,Float64}, size(df_gene, 1))
    for i = 1:size(df_gene, 1)
        source = df_gene.Source[i]
        target = df_gene.Target[i]
        if source ∉ keys(location_dict) || target ∉ keys(location_dict)
            distance[i] = missing
        elseif location_dict[source][1] != location_dict[target][1]
            distance[i] = -1
        else
            source_start = location_dict[source][3]
            target_start = location_dict[target][3]
            distance[i] = abs(source_start - target_start)
        end
    end
    df_gene[!, :Distance] = distance
    df_gene[!, :same_chromosome] = distance .>= 0
    df_gene[!, :link_type] .= "gene-gene"
    # prepare methylation dataframe
    df_methylation[!, :same_chromosome] .= true
    # return integration of gene and methylation network
    return append!(df_gene, df_methylation)
end


function metagraph_from_dataframe(
    graph_type,
    df::DataFrame,
    origin::Symbol,
    destination::Symbol,
    weight::Symbol = Symbol(),
    edge_attributes::Union{Vector{Symbol},Symbol} = Vector{Symbol}(),
    vertex_attributes::DataFrame = DataFrame(),
    vertex_id_col::Symbol = Symbol(),
)

    # Map node names to vertex IDs
    nodes = sort!(unique!([df[:, origin]; df[:, destination]]))
    vertex_names = DataFrame(name = nodes, vertex_id = eachindex(nodes))

    # Merge in to original
    for c in [origin, destination]
        temp = rename(vertex_names, :vertex_id => Symbol(c, :_id))
        df = innerjoin(df, temp; on = c => :name)
    end

    # Merge additional attributes to names
    if vertex_attributes != DataFrame()
        idsym =
            vertex_id_col == Symbol() ? first(propertynames(vertex_attributes)) :
            vertex_id_col
        vertex_names = leftjoin(vertex_names, vertex_attributes, on = :name => idsym)
    end

    # Create Graph
    mg = graph_type(nrow(vertex_names))
    for r in eachrow(df)
        add_edge!(mg, r[Symbol(origin, :_id)], r[Symbol(destination, :_id)])
    end

    # Set vertex names and attributes
    attr_names = propertynames(vertex_names[!, Not(:vertex_id)])
    for r in eachrow(vertex_names)
        set_props!(mg, r[:vertex_id], Dict([a => r[a] for a in attr_names]))
    end

    # Set edge attributes
    if typeof(edge_attributes) == Symbol
        edge_attributes = Vector{Symbol}([edge_attributes])
    end

    origin_id = Symbol(origin, :_id)
    destination_id = Symbol(destination, :_id)

    for e in edge_attributes
        for r in eachrow(df)
            set_prop!(mg, Edge(r[origin_id], r[destination_id]), e, r[e])
        end
    end

    # Weight
    if weight != Symbol()
        for r in eachrow(df)
            set_prop!(mg, Edge(r[origin_id], r[destination_id]), :weight, r[weight])
        end
    end

    # Set name as index (Issue #9)
    set_indexing_prop!(mg, :name)

    return mg
end

function find_community(g::AbstractMetaGraph, γ = 0.1)
    link_df = DataFrame(from = Int[], to = Int[], weight = Float64[])
    node_df = DataFrame(id = Int[], label = String[], importance = Float64[])
    for v in vertices(g)
        push!(node_df, (v, get_prop(g, v, :name), 1.0))
    end
    for e in edges(g)
        push!(link_df, (e.src, e.dst, 1.0))
    end
    jc = juliac.JuliaCommunityInstance(
        link_df,
        nodes = node_df,
        node_label_field = "label",
        edge_weighted = false,
        node_weighted = false,
        is_directed = false,
        to_summarise_graph = false
    )
    jc.γ = γ
    juliac.discover_communities(jc)
    return jc
end


function create_jc_instance(g::AbstractMetaGraph)
    link_df = DataFrame(from = Int[], to = Int[], weight = Float64[])
    node_df = DataFrame(id = Int[], label = String[], importance = Float64[])
    for v in vertices(g)
        push!(node_df, (v, get_prop(g, v, :name), 1.0))
    end
    for e in edges(g)
        push!(link_df, (e.src, e.dst, 1.0))
    end
    jc = juliac.JuliaCommunityInstance(
        link_df,
        nodes = node_df,
        node_label_field = "label",
        edge_weighted = false,
        node_weighted = false,
        is_directed = false,
        to_summarise_graph = false
    )
    return jc
end

function metagraph2nodeset_df(g)
    df = DataFrame()
    for v in vertices(g)
        push!(df, props(g, v), cols = :union)
    end
    rename!(df, Dict(:name => :id))
    return df
end

function metagraph2edge_df(g)
    df = DataFrame()
    for e in edges(g)
        input_dict = props(g, e)
        input_dict[:Source] = props(g, src(e))[:name]
        input_dict[:Target] = props(g, dst(e))[:name]

        push!(df, input_dict, cols = :union)
    end
    return df
end
##
# Select the edges in top 1% of the PIDC rank
df_gene_tumor = get_df_network(tumor_input_network_path, threshold = 0.99)

df_gene_loc = CSV.read(input_gene_loc_path, DataFrame, select = ["ensembl_gene_id", "external_gene_name", "seqnames", "start", "end"])
df_methylation = get_df_methylation(hypermethylation_path = input_hypermethylation_path, hypomethylation_path = input_hypomethylation_path)
gene_node_set = Set([df_gene_tumor.Source; df_gene_tumor.Target; df_methylation.Target])
gene_location_dict = get_gene_location_dictionary(gene_node_set, df_gene_loc)

# integrate the DMCs network to DEGs network
add_methylation_node!(gene_location_dict, df_methylation)
df_node_set = create_node_set_dataframe(gene_location_dict)
df_network = create_gene_methylation_network_dataframe(df_gene_tumor, df_methylation, gene_location_dict)


##

mg_ = metagraph_from_dataframe(MetaGraph, df_network, :Source, :Target, Symbol(), [:Distance, :same_chromosome, :link_type], df_node_set, :id)

# Extract the giant component of the graph
greatest_component_vertices = sort(connected_components(mg_), by = length, rev = true)[1]
mg, _ = induced_subgraph(mg_, greatest_component_vertices)


# Run this code to determine the optimal resolution value
#jc = create_jc_instance(mg)
#juliac.optimise_resolution(jc,γ_from=0.0001, γ_end=0.05, γ_step=0.0002)

# Community detection using leiden_community.jl
jc = find_community(mg, 0.0085)

# Add leiden membership to node properties
for v in vertices(mg)
    set_prop!(mg, v, :leiden_membership, jc.membership_vector[v])
end


##
# Calculate centrality
between_cent = Parallel.betweenness_centrality(mg)
k_core_cent = core_number(mg)
degree_cent = degree(mg)

# Add centrality to node properties
for v in vertices(mg)
    set_props!(mg, v, Dict(:betweenness => between_cent[v],
        :k_core => k_core_cent[v],
        :degree => degree_cent[v]))
end



##
# Convert network object to nodes dataframe and edges dataframe
new_df_node_set = metagraph2nodeset_df(mg)
new_df_network = metagraph2edge_df(mg)
new_df_network[!, :Target_member] .= 0
new_df_network[!, :Source_member] .= 0
##
# Add properti of leiden membership for source node and target node in edges dataframe.
# Source and target do not mean the edges have direction (because the network is undirected).
# It only serve to simplyfy the naming of columns.
target_source_member = Vector{Set{Int64}}(undef, size(new_df_network, 1))
for i = 1:size(new_df_network, 1)
    target = new_df_network[i, :Target]
    source = new_df_network[i, :Source]
    target_member = filter(:id => x -> x == target, new_df_node_set)[!, :leiden_membership][1]
    source_member = filter(:id => x -> x == source, new_df_node_set)[!, :leiden_membership][1]
    target_source_member[i] = Set([target_member, source_member])
    new_df_network[i, :Target_member] = target_member
    new_df_network[i, :Source_member] = source_member
end
##
#Save the result of node dataframe and edges dataframe
CSV.write(tumor_output_network_path, new_df_network)
CSV.write(tumor_output_gene_loc_path, new_df_node_set)
##
