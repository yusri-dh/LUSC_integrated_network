using Distributed
addprocs(20,exeflags=`--project=$(Base.active_project())`)
@everywhere using NetworkInference

DEGs_filepath = "./TCGA-LUSC_Tumor.tsv"

function create_network(input)
    dir_name = dirname(input)
    file_name = basename(input)
    nodes = get_nodes(input)
    inferred_network = InferredNetwork(PIDCNetworkInference(), nodes)
    output = joinpath(dir_name, "output_" * file_name)
    println(output)
    write_network_file(output, inferred_network)
    println("OK")
    return inferred_network
end

inferred_network = create_network(DEGs_filepath)
