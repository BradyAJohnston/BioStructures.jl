export StructureArray, ChainArray, ResidueArray, AtomArrray, find_offsets, getvalue, chains, residues, resdiues, resnames, coords, resnameselector

struct StructureArray
    n_atoms::UInt64
    chain_offsets::Vector{UInt64}
    residue_offsets::Vector{UInt64}
    data::Dict
    coord::Matrix{Float64}

    function StructureArray(data::Dict)
        coord = Matrix(transpose(hcat(data["Cartn_x"], data["Cartn_y"], data["Cartn_z"])))
        [pop!(data, key) for key in ["Cartn_x", "Cartn_y", "Cartn_z"]]
        new(
            length(data["auth_asym_id"]),
            find_offsets(data["auth_asym_id"]),
            find_offsets(data["auth_seq_id"]),
            data,
            coord
        )
    end
end

Base.show(io::IO, struc::StructureArray) = print(io, "StructureArray with ", struc.n_atoms, " atoms, ", length(struc.chain_offsets), " chains, ", length(struc.residue_offsets), " residues")

struct ChainArray
    struc_array::StructureArray
    start_idx::UInt64
    end_idx::UInt64
    chain_id::String
    residue_offsets::Vector{UInt64}
end

Base.show(io::IO, chain::ChainArray) = print(io, "ChainArray(chain ", chain.chain_id, ", ", length(chain.residue_offsets), " residues, indices ", chain.start_idx, ":", chain.end_idx, ")")

struct ResidueArray
    struc_array::StructureArray
    start_idx::UInt64
    end_idx::UInt64
    chain_id::String
    res_id::String

    function ResidueArray(struc_array::StructureArray, start_idx::UInt64, end_idx::UInt64, chain_id::String, res_id::String)
        new(struc_array, start_idx, end_idx, chain_id, res_id)
    end
end

Base.show(io::IO, residue::ResidueArray) = print(io, "ResidueArray(", residue.res_id, " in chain ", residue.chain_id, ", ", residue.end_idx - residue.start_idx + 1, " atoms, indices ", residue.start_idx, ":", residue.end_idx, ")")

struct AtomArray
    struc_array::StructureArray
    idx::UInt64

    function AtomArrray(struc_array::StructureArray, idx::UInt64)
        new(struc_array, idx)
    end
end

function find_offsets(arr::Union{Vector{Int32},Vector{String}})
    offsets = findall([true, (arr[1:end-1] .!= arr[2:end])...])
    return sort(offsets)
end

function get_chain(struc::StructureArray, chain_id::String)
    for (i, offset) in enumerate(struc.chain_offsets)
        if struc.data["auth_asym_id"][offset] == chain_id
            if i == length(struc.chain_offsets)
                end_idx = struc.n_atoms
            else
                end_idx = struc.chain_offsets[i+1] - 1
            end
            residue_offset_mask = (struc.residue_offsets .>= offset) .& (struc.residue_offsets .<= end_idx)
            return ChainArray(struc, offset, end_idx, chain_id, struc.residue_offsets[residue_offset_mask])
        end
    end
    throw(ArgumentError("Chain ID $chain_id not found"))
end

function get_chain(struc::StructureArray, index::Int)
    if index < 1 || index > length(struc.chain_offsets)
        throw(ArgumentError("Chain index out of bounds"))
    end
    offset = struc.chain_offsets[index]
    if index == length(struc.chain_offsets)
        end_idx = struc.n_atoms
    else
        end_idx = struc.chain_offsets[index+1] - 1
    end
    residue_offset_mask = (struc.residue_offsets .>= offset) .& (struc.residue_offsets .<= end_idx)
    chain_id = struc.data["auth_asym_id"][offset]
    return ChainArray(struc, offset, end_idx, chain_id, struc.residue_offsets[residue_offset_mask])
end

Base.getindex(struc::StructureArray, i::Int) = get_chain(struc, i)
Base.getindex(struc::StructureArray, chain_id::String) = get_chain(struc, chain_id)

function chains(struc::StructureArray)
    return [get_chain(struc, i) for i in 1:length(struc.chain_offsets)]
end

function residues(struc::StructureArray)
    return [get_residue(struc, i) for i in 1:length(struc.residue_offsets)]
end

function get_residue(chain::ChainArray, res_num::Int)
    if res_num < 1 || res_num > length(chain.residue_offsets) - 1
        throw(ArgumentError("Residue number out of bounds"))
    end
    start_idx = chain.residue_offsets[res_num]
    end_idx = chain.residue_offsets[res_num+1] - 1
    return ResidueArray(chain.struc_array, start_idx, end_idx, chain.chain_id, chain.struc_array.data["label_comp_id"][start_idx])
end

function get_residue(struc::StructureArray, index::Int)
    if index < 1 || index > length(struc.residue_offsets)
        throw(ArgumentError("Residue index out of bounds"))
    end
    start_idx = struc.residue_offsets[index]
    if index == length(struc.residue_offsets)
        end_idx = UInt(length(struc.data["auth_asym_id"]))
    else
        end_idx = struc.residue_offsets[index+1] - 1
    end
    chain_id = struc.data["auth_asym_id"][start_idx]
    return ResidueArray(struc, start_idx, end_idx, chain_id, struc.data["label_comp_id"][start_idx])
end

Base.getindex(chain::ChainArray, i::Int) = get_residue(chain, i)

function residues(chain::ChainArray)
    residues = []
    offsets = chain.residue_offsets
    for (i, offset) in enumerate(offsets[1:end-1])
        push!(residues, ResidueArray(chain.struc_array, offset, offsets[i+1] - 1, chain.chain_id, chain.struc_array.data["label_comp_id"][offset]))
    end
    return residues
end


function get_residue(chain::ChainArray, res_num::Int)
    offsets = offsets(chain)
    if res_num < 1 || res_num > length(offsets) - 1
        throw(ArgumentError("Residue number out of bounds"))
    end
    start_idx = offsets[res_num]
    end_idx = offsets[res_num+1] - 1
    return ResidueArray(chain.struc_array, start_idx, end_idx, chain.chain_id, chain.struc_array.data["label_comp_id"][start_idx])
end

function getvalue(chain::ChainArray, value::String)
    return chain.struc_array.data[value][chain.start_idx:chain.end_idx]
end

function getvalue(struc::StructureArray, value::String)
    return struc.data[value]
end

function getvalue(residue::ResidueArray, value::String)
    return residue.struc_array.data[value][residue.start_idx:residue.end_idx]
end

# Base.iterate(chain::ChainArray, i) = iterate(get_residue(chain, i), i)

Base.iterate(struc::StructureArray, i) = iterate(get_chain(struc, i), i)
Base.iterate(struc::StructureArray) = iterate([get_chain(struc, i) for i in 1:length(struc.chain_offsets)], 1)
Base.iterate(chain::ChainArray, i) = iterate(residues(chain), i)
Base.length(struc::StructureArray) = length(struc.chain_offsets)
Base.length(chain::ChainArray) = length(chain.start_idx)
Base.length(residue::ResidueArray) = length(residue.start_idx)

function coords(struc::StructureArray)

    return struc.coord
end
function coords(chain::ChainArray)
    return chain.struc_array.coord[:, chain.start_idx:chain.end_idx]
end
function coords(residue::ResidueArray)
    return residue.struc_array.coord[:, residue.start_idx:residue.end_idx]
end
function coords(atom::AtomArray)
    return atom.struc_array.coord[:, atom.idx]
end


function resnames(struc::StructureArray)
    return struc.data["label_comp_id"][struc.residue_offsets]
end
function resnames(chain::ChainArray)
    return chain.struc_array.data["label_comp_id"][chain.residue_offsets]
end

function resnameselector(el::StructureArray, resnames::String)
    return resnames(el) .==  resname
end

function resnameselector(el::StructureArray, resnames::Vector{String})
    return reduce(+, [resnames(el) .== r for r in resnames])
end

function chainids(struc::StructureArray)
    return struc.data["auth_asym_id"][struc.chain_offsets]
end
