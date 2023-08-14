module XCrySDenStructureFormat

using AtomsBase
using Unitful
using UnitfulAtomic
using PeriodicTable: PeriodicTable
using StaticArrays

const LENGTH_UNIT = u"Å"
const FORCE_UNIT = u"Eh_au" / u"Å"

export load_xsf
export save_xsf
include("fileio.jl")

end # module XSD
