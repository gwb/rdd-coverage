import LibGEOS
import GeoJSON
import Base.convert
using DataFrames
import CSV
import JSON

import LibGEOS: LineString, MultiLineString
import GeoInterface

# type conversion
const BorderType = Union{MultiLineString, LineString}
BorderType(ls::GeoInterface.LineString) = LineString(ls)
BorderType(ls::GeoInterface.MultiLineString) = MultiLineString(ls)

const SchDistr = Int
const RegionType = Union{GeoInterface.AbstractPolygon,
                           GeoInterface.AbstractMultiPolygon}
convert(::Type{RegionType}, p::GeoInterface.Polygon) = LibGEOS.Polygon(p)
convert(::Type{RegionType}, p::GeoInterface.MultiPolygon) = LibGEOS.MultiPolygon(p)

const SALE_PRICE = Symbol("SALE PRICE")
const SQFT = Symbol("GROSS SQUARE FEET")
const BUILDING_CLASS_AT_TIME_OF_SALE = Symbol("BUILDING CLASS AT TIME OF SALE")
const BUILDING_CLASS_CATEGORY = Symbol("BUILDING CLASS CATEGORY")
const TAX_CLASS_AT_TIME_OF_SALE = Symbol("TAX CLASS AT TIME OF SALE")

function read_borders()
    SchoolDistrict_borders=GeoJSON.parsefile(
            "NYC_data/processed"
            *"/SchoolDistrict_borders/SchoolDistrict_borders.json")
    borders_dict = Dict{Tuple{SchDistr,SchDistr}, 
                        Union{GeoInterface.AbstractMultiLineString,
                              GeoInterface.AbstractLineString}}()

    for f in SchoolDistrict_borders.features
        schdist1 = f.properties["SchoolDistrict1"]
        schdist2 = f.properties["SchoolDistrict2"]
        border = f.geometry
        borders_dict[schdist1,schdist2] = BorderType(border)
    end
    return borders_dict
end

function read_distr_shapes()
    json=GeoJSON.parsefile("NYC_data/processed"
                           *"/SchoolDistrict_shapes/SchoolDistrict_shapes.json")
    schdistr_shape_dict = Dict{SchDistr, RegionType}()
    for feature in json.features
        schdistr = convert(SchDistr, feature.properties["SchoolDist"])
        shape = feature.geometry
        schdistr_shape_dict[schdistr] = shape
    end
    return schdistr_shape_dict
end

function read_distr_reprpoints()
    json = JSON.parsefile("NYC_data/processed"
                          *"/SchoolDistrict_reprpoints/schdistr_representative.json")
    schdistr_repr_dict = Dict{SchDistr, Tuple{Float64, Float64}}()
    for schdistr_int in keys(json)
        schdistr = parse(SchDistr, schdistr_int)
        value = json[schdistr_int]
        v1 = value[1]::Float64
        v2 = value[2]::Float64
        schdistr_repr_dict[schdistr] = (v1, v2)
    end
    return schdistr_repr_dict
end

function read_processed_sales()
    NYC_sales=CSV.read("NYC_data/processed/NYC_sales.csv", DataFrame,
                   types=Dict("TAX CLASS AT PRESENT" => Union{Missings.Missing, String},
                              "TAX CLASS AT TIME OF SALE" => Union{Missings.Missing, String}),
                   weakrefstrings=false,
                   nullable=true,
                   )
    nyc_schdistrs = NYC_sales[:SchDistr]
    schd_strings = [ismissing(sd)?missing:dec(sd,2) for sd in  nyc_schdistrs]
    str_schdistrs = CategoricalVector(schd_strings)
    NYC_sales[:SchDistr] = str_schdistrs
    # categorical variables
    categorical!(NYC_sales, Symbol["BOROUGH",
        "BUILDING CLASS CATEGORY",
        "BUILDING CLASS AT TIME OF SALE",
        "ZIP CODE",
        "TAX CLASS AT PRESENT",
        "TAX CLASS AT TIME OF SALE",
        "NEIGHBORHOOD",
        ])
    sort!(NYC_sales; cols=:SchDistr)
    return NYC_sales
end

const DWELLINGS_DICT = Dict(
    "01  ONE FAMILY DWELLINGS"=>true,
    "02  TWO FAMILY DWELLINGS"=>true,
    "03  THREE FAMILY DWELLINGS"=>true,
    "04  TAX CLASS 1 CONDOS"=>false,
    "05  TAX CLASS 1 VACANT LAND"=>false,
    "06  TAX CLASS 1 - OTHER"=>false,
    "07  RENTALS - WALKUP APARTMENTS"=>false,
    "08  RENTALS - ELEVATOR APARTMENTS"=>false,
    "09  COOPS - WALKUP APARTMENTS"=>false,
    "10  COOPS - ELEVATOR APARTMENTS"=>false,
    "11  SPECIAL CONDO BILLING LOTS"=>false,
    "11A CONDO-RENTALS"=>false,
    "12  CONDOS - WALKUP APARTMENTS"=>false, # why are these duplicated?
    "13  CONDOS - ELEVATOR APARTMENTS"=>false,
    "14  RENTALS - 4-10 UNIT"=>false,
    "15  CONDOS - 2-10 UNIT RESIDENTIAL"=>false,
    "16  CONDOS - 2-10 UNIT WITH COMMERCIAL UNIT"=>false,
    "17  CONDO COOPS"=>false,
    "18  TAX CLASS 3 - UNTILITY PROPERTIES"=>false,
    "21  OFFICE BUILDINGS"=>false,
    "22  STORE BUILDINGS"=>false,
    "23  LOFT BUILDINGS"=>false,
    "25  LUXURY HOTELS"=>false,
    "26  OTHER HOTELS"=>false,
    "27  FACTORIES"=>false,
    "28  COMMERCIAL CONDOS"=>false,
    "29  COMMERCIAL GARAGES"=>false,
    "30  WAREHOUSES"=>false,
    "31  COMMERCIAL VACANT LAND"=>false,
    "32  HOSPITAL AND HEALTH FACILITIES"=>false,
    "33  EDUCATIONAL FACILITIES"=>false,
    "34  THEATRES"=>false,
    "35  INDOOR PUBLIC AND CULTURAL FACILITIES"=>false,
    "36  OUTDOOR RECREATIONAL FACILITIES"=>false,
    "37  RELIGIOUS FACILITIES"=>false,
    "38  ASYLUMS AND HOMES"=>false,
    "41  TAX CLASS 4 - OTHER"=>false,
    "42  CONDO CULTURAL/MEDICAL/EDUCATIONAL/ETC"=>false,
    "43  CONDO OFFICE BUILDINGS"=>false,
    "44  CONDO PARKING"=>false,
    "45  CONDO HOTELS"=>false,
    "46  CONDO STORE BUILDINGS"=>false,
    "47  CONDO NON-BUSINESS STORAGE"=>false,
    "48  CONDO TERRACES/GARDENS/CABANAS"=>false,
    "49  CONDO WAREHOUSES/FACTORY/INDUS"=>false,
    )
const RESIDENTIAL_DICT = Dict(
    "01  ONE FAMILY DWELLINGS"=>true,
    "02  TWO FAMILY DWELLINGS"=>true,
    "03  THREE FAMILY DWELLINGS"=>true,
    "04  TAX CLASS 1 CONDOS"=>true,
    "05  TAX CLASS 1 VACANT LAND"=>false,
    "06  TAX CLASS 1 - OTHER"=>false,
    "07  RENTALS - WALKUP APARTMENTS"=>false,
    "08  RENTALS - ELEVATOR APARTMENTS"=>false,
    "09  COOPS - WALKUP APARTMENTS"=>true,
    "10  COOPS - ELEVATOR APARTMENTS"=>true,
    "11  SPECIAL CONDO BILLING LOTS"=>false,
    "11A CONDO-RENTALS"=>false,
    "12  CONDOS - WALKUP APARTMENTS"=>true, # why are these duplicated?
    "13  CONDOS - ELEVATOR APARTMENTS"=>true,
    "14  RENTALS - 4-10 UNIT"=>false,
    "15  CONDOS - 2-10 UNIT RESIDENTIAL"=>true,
    "16  CONDOS - 2-10 UNIT WITH COMMERCIAL UNIT"=>false,
    "17  CONDO COOPS"=>true,
    "18  TAX CLASS 3 - UNTILITY PROPERTIES"=>false,
    "21  OFFICE BUILDINGS"=>false,
    "22  STORE BUILDINGS"=>false,
    "23  LOFT BUILDINGS"=>false,
    "25  LUXURY HOTELS"=>false,
    "26  OTHER HOTELS"=>false,
    "27  FACTORIES"=>false,
    "28  COMMERCIAL CONDOS"=>false,
    "29  COMMERCIAL GARAGES"=>false,
    "30  WAREHOUSES"=>false,
    "31  COMMERCIAL VACANT LAND"=>false,
    "32  HOSPITAL AND HEALTH FACILITIES"=>false,
    "33  EDUCATIONAL FACILITIES"=>false,
    "34  THEATRES"=>false,
    "35  INDOOR PUBLIC AND CULTURAL FACILITIES"=>false,
    "36  OUTDOOR RECREATIONAL FACILITIES"=>false,
    "37  RELIGIOUS FACILITIES"=>false,
    "38  ASYLUMS AND HOMES"=>false,
    "41  TAX CLASS 4 - OTHER"=>false,
    "42  CONDO CULTURAL/MEDICAL/EDUCATIONAL/ETC"=>false,
    "43  CONDO OFFICE BUILDINGS"=>false,
    "44  CONDO PARKING"=>false,
    "45  CONDO HOTELS"=>false,
    "46  CONDO STORE BUILDINGS"=>false,
    "47  CONDO NON-BUSINESS STORAGE"=>false,
    "48  CONDO TERRACES/GARDENS/CABANAS"=>false,
    "49  CONDO WAREHOUSES/FACTORY/INDUS"=>false,
    )

function filter_sales(NYC_sales::DataFrame)

    NYC_sales[:logSalePricePerSQFT] = map(log, NYC_sales[SALE_PRICE]) .- map(log, NYC_sales[SQFT])
    believable = zeros(Bool, size(NYC_sales,1))
    removed=zeros(Int, 12)
    for i in 1:size(NYC_sales,1)
        if ismissing(NYC_sales[i,SALE_PRICE])
            # remove sales without a sale price
            removed[1] += 1
        elseif ismissing(NYC_sales[i,:SchDistr])
            # remove properties without a school district
            removed[2] += 1
        elseif ismissing(NYC_sales[i,BUILDING_CLASS_AT_TIME_OF_SALE])
            # remove sales with missing covariates
            removed[3] += 1
        elseif ismissing(NYC_sales[i,BUILDING_CLASS_CATEGORY])
            # remove sales with missing covariates
            removed[4] += 1
        elseif !DWELLINGS_DICT[NYC_sales[i,BUILDING_CLASS_CATEGORY]]
            # remove non-residential properties
            # in fact, remove things that aren't dwellings (houses?)
            removed[5] += 1
        elseif ismissing(NYC_sales[i,BUILDING_CLASS_CATEGORY])
            # remove sales with missing covariates
            removed[6] += 1
        elseif ismissing(NYC_sales[i,TAX_CLASS_AT_TIME_OF_SALE])
            # remove sales with missing covariates
            removed[7] += 1
        elseif ismissing(NYC_sales[i,SQFT])
            # remove sales with missing GROSS SQUARE FEET information
            removed[8] += 1
        elseif NYC_sales[i,SQFT]<100.0
            # remove properties smaller to 100sqft (seems unlikely to be real)
            removed[9] += 1
        elseif NYC_sales[i,:logSalePricePerSQFT] < 3.0
            # that's too cheap (remove outliers)
            removed[10] += 1
        elseif NYC_sales[i,:logSalePricePerSQFT] > 8.0
            # that's too expensive (remove outliers)
            removed[11] += 1
        elseif ismissing(NYC_sales[i,:XCoord])
            # remove properties with failed geocoding
            removed[12] += 1
        elseif ismissing(NYC_sales[i,:YCoord])
            # remove properties with failed geocoding
            removed[13] += 1
        else
            # otherwise keep
            believable[i] = true
        end
    end

    filtered = copy(NYC_sales[believable,:])
    return Dict(
        :filtered => filtered,
        :believable => believable,
        :removed => removed,
        )
end

function read_sentinels()
    sentinels_json=GeoJSON.parsefile("NYC_data/processed/SchoolDistrict_borders/SchoolDistrict_sentinels.json")
    sentinels=Dict{Tuple{String,String},GeoInterface.MultiPoint}()
    for f in sentinels_json.features
        key = (dec(f.properties["SchoolDistrict1"],2), dec(f.properties["SchoolDistrict2"],2))
        sentinels[key] = f.geometry
    end
    return sentinels
end

function sales_dicts(NYC_sales::DataFrame)
    schdistrs = sort(NYC_sales[:SchDistr].pool.levels)
    schdistr_indices = Dict{String,Vector{Int}}()
    schdistrs_col = NYC_sales[:SchDistr]
    for name in schdistrs
        indices = find(schdistrs_col.refs .== find(schdistrs_col.pool.index .== name)[1])
        nobsv_schdistr = length(indices)
        @printf("District %s has %d sales\n", name, nobsv_schdistr)
        schdistr_indices[name] = indices
    end
    Y_dict=Dict{String, Vector{Float64}}()
    X_dict=Dict{String, Array{Float64,2}}()
    for name in schdistrs
        Y_dict[name] = NYC_sales[schdistr_indices[name], :logSalePricePerSQFT]
        X_dict[name] = NYC_sales[schdistr_indices[name],[:XCoord, :YCoord]]
    end
    return Dict(
            :schdistr_indices => schdistr_indices,
            :schdistrs => schdistrs,
            :X_dict => X_dict,
            :Y_dict => Y_dict
            )
end


