using Enzyme

struct Edge
    from::Integer
    to::Integer
end

struct Triangle
    p1::Integer
    p2::Integer
    p3::Integer
end

struct Surface
    edges::Vector{Edge}
    triangles::Vector{Triangle}
end


# X, Y grid positions for points
positions = Dict{Integer, Vector{Real}}(
    1=>[-0.5, -0.5],
    2=>[0.5, -0.5],
    3=>[0.5, 0.5],
    4=>[-0.5, 0.5],
    5=>[-0.3, -0.5],
    6=>[-0.1, -0.5],
    7=>[0.1, -0.5],
    8=>[0.3, -0.5],
    9=>[0.5, -0.3],
    10=>[0.5, -0.1],
    11=>[0.5, 0.1],
    12=>[0.5, 0.3],
    13=>[0.3, 0.5],
    14=>[0.1, 0.5],
    15=>[-0.1, 0.5],
    16=>[-0.3, 0.5],
    17=>[-0.5, 0.3],
    18=>[-0.5, 0.1],
    19=>[-0.5, -0.1],
    20=>[-0.5, -0.3],
    21=>[-0.3, -0.3],
    22=>[-0.1, -0.3],
    23=>[0.1, -0.3],
    24=>[0.3, -0.3],
    25=>[0.3, -0.1],
    26=>[0.3, 0.1],
    27=>[0.3, 0.3],
    28=>[0.1, 0.3],
    29=>[-0.1, 0.3],
    30=>[-0.3, 0.3],
    31=>[-0.3, 0.1],
    32=>[-0.3, -0.1],
    33=>[-0.1, -0.1],
    34=>[0.1, -0.1],
    35=>[0.1, 0.1],
    36=>[-0.1, 0.1],
)

# Height Function
height(x, y) = 500 - (x^2 + y^2) * 10

# Creating points on the surface
xpoints = Float64[]
ypoints = Float64[]

zpoints = Float64[]
for i in eachindex(1:length(positions))
    push!(xpoints, positions[i][1])
    push!(ypoints, positions[i][2])
    push!(zpoints, height(positions[i][1], positions[i][2]))
end

# Setup Edge connectivity
edges = Edge[
    Edge(1, 5),         #1
    Edge(5, 6),         #2
    Edge(6, 7),         #3
    Edge(7, 8),         #4
    Edge(8, 2),         #5
    Edge(2, 9),         #6
    Edge(9, 10),         #7
    Edge(10, 11),         #8
    Edge(11, 12),         #9
    Edge(12, 3),         #10
    Edge(3, 13),         #11
    Edge(13, 14),         #12
    Edge(14, 15),         #13
    Edge(15, 16),         #14
    Edge(16, 4),         #15
    Edge(4, 17),         #16
    Edge(17, 18),         #17
    Edge(18, 19),         #18
    Edge(19, 20),         #19
    Edge(20, 1),         #20
    Edge(20, 21),         #21
    Edge(21, 22),         #22
    Edge(22, 23),         #23
    Edge(23, 24),         #24
    Edge(24, 9),         #25
    Edge(19, 32),         #26
    Edge(32, 33),         #27
    Edge(33, 34),         #28
    Edge(34, 25),         #29
    Edge(25, 10),         #30
    Edge(18, 31),         #31
    Edge(31, 36),         #32
    Edge(36, 35),         #33
    Edge(35, 26),         #34
    Edge(26, 11),         #35
    Edge(17, 30),         #36
    Edge(30, 29),         #37
    Edge(29, 28),         #38
    Edge(28, 27),         #39
    Edge(27, 12),         #40
    Edge(5, 21),         #41
    Edge(21, 32),         #42
    Edge(32, 31),         #43
    Edge(31, 30),         #44
    Edge(30, 16),         #45
    Edge(6, 22),         #46
    Edge(22, 33),         #47
    Edge(33, 36),         #48
    Edge(36, 29),         #49
    Edge(29, 15),         #50
    Edge(7, 23),         #51
    Edge(23, 34),         #52
    Edge(34, 35),         #53
    Edge(35, 28),         #54
    Edge(28, 14),         #55
    Edge(8, 24),         #56
    Edge(24, 25),         #57
    Edge(25, 26),         #58
    Edge(26, 27),         #59
    Edge(27, 13),         #60
    Edge(1, 21),         #61
    Edge(5, 22),         #62
    Edge(6, 23),         #63
    Edge(7, 24),         #64
    Edge(8, 9),         #65
    Edge(20, 32),         #66
    Edge(21, 33),         #67
    Edge(22, 34),         #68
    Edge(23, 25),         #69
    Edge(24, 10),         #70
    Edge(19, 31),         #71
    Edge(32, 36),         #72
    Edge(33, 35),         #73
    Edge(34, 26),         #74
    Edge(25, 11),         #75
    Edge(18, 30),         #76
    Edge(31, 29),         #77
    Edge(36, 28),         #78
    Edge(35, 27),         #79
    Edge(26, 12),         #80
    Edge(17, 16),         #81
    Edge(30, 15),         #82
    Edge(29, 14),         #83
    Edge(28, 13),         #84
    Edge(27, 3),         #85
]

# Setup Face connectivity
triangles = Triangle[
    Triangle(1,5,21),
    Triangle(5,6,22),
    Triangle(6,7,23),
    Triangle(7,8,24),
    Triangle(8,2,9),
    Triangle(1,21,20),
    Triangle(5,22,21),
    Triangle(6,23,22),
    Triangle(7,24,23),
    Triangle(8,9,24),
    Triangle(20,21,32),
    Triangle(21,22,33),
    Triangle(22,23,34),
    Triangle(23,24,25),
    Triangle(24,9,10),
    Triangle(20,32,19),
    Triangle(21,33,32),
    Triangle(22,34,33),
    Triangle(23,25,34),
    Triangle(24,10,9),
    Triangle(19,32,31),
    Triangle(32,33,36),
    Triangle(33,34,35),
    Triangle(35,25,26),
    Triangle(25,10,11),
    Triangle(19,31,18),
    Triangle(32,36,31),
    Triangle(33,35,36),
    Triangle(34,26,35),
    Triangle(25,11,26),
    Triangle(18,31,30),
    Triangle(31,36,29),
    Triangle(36,35,28),
    Triangle(35,26,27),
    Triangle(26,11,12),
    Triangle(18,30,17),
    Triangle(31,29,30),
    Triangle(36,28,29),
    Triangle(35,27,28),
    Triangle(26,12,27),
    Triangle(17,30,16),
    Triangle(30,29,15),
    Triangle(29,28,14),
    Triangle(28,27,13),
    Triangle(27,12,3),
    Triangle(17,16,4),
    Triangle(30,15,16),
    Triangle(29,14,15),
    Triangle(28,13,14),
    Triangle(27,3,13),
]


# Create surface
s = Surface(edges, triangles)


function calculate_area(s::Surface, xpoints, ypoints, zpoints)
    area = zero(eltype(xpoints))
    for triangle in s.triangles
        x1 = xpoints[triangle.p1]
        y1 = ypoints[triangle.p1]
        z1 = zpoints[triangle.p1]
        x2 = xpoints[triangle.p2]
        y2 = ypoints[triangle.p2]
        z2 = zpoints[triangle.p2]
        x3 = xpoints[triangle.p3]
        y3 = ypoints[triangle.p3]
        z3 = zpoints[triangle.p3]
        u1 = x2 - x1
        u2 = y2 - y1
        u3 = z2 - z1
        v1 = x3 - x1
        v2 = y3 - y1
        v3 = z3 - z1
        area = area + 0.5 * sqrt((u2*v3 - u3*v2)^2 + (u3*v1-u1*v3)^2 + (u1*v2-u2*v1)^2)
    end
    return area
end

@show calculate_area(s, xpoints, ypoints, zpoints)

#@show Enzyme.gradient(Reverse, calculate_area, Const(s), Const(xpoints), Const(ypoints), zpoints)

using Zygote, BenchmarkTools

@btime Zygote.gradient(calculate_area, s, xpoints, ypoints, zpoints)

#@btime Enzyme.gradient(Reverse, $calculate_area, $(Const(s)), $(Const(xpoints)), $(Const(ypoints)), $(zpoints))
#  62.000 μs (361 allocations: 41.61 KiB)
# (nothing, nothing, nothing, [-0.1408721450121023, -0.14087214501209985, -0.1408721450121023, -0.14087214501209985, -0.13333333333332764, -0.143672232115815, -0.1436722321158125, -0.13333333333333264, -0.2666666666666653, -0.010338898782480455  …  0.09049890913642407, 0.13554687099016932, 0.13554687099016688, 0.09049890913643148, 0.13554687099016688, 0.13554687099016935, 0.3332906247936249, 0.2362763747790916, 0.4303048748081582, 0.3332906247936249])

@btime Zygote.gradient($calculate_area, $s, $xpoints, $ypoints, $zpoints)
#  1.043 ms (17118 allocations: 635.86 KiB)
# (nothing, [-0.008804509063236488, 0.008804509063256116, 0.008804509063236154, -0.008804509063256116, 0.7691150024933462, 0.2715520029479166, -0.2715520029479066, -0.7691150024933762, 0.8024483358267235, 0.7235145845699962  …  0.5419365898655986, 0.27741340607116555, -0.27741340607115533, -0.5419365898656485, -0.7529730203651925, -0.7529730203652125, -0.4706492344625691, 0.4949027969662021, 0.44639567195893615, -0.4706492344625691], [-0.008804509063236488, -0.008804509063256116, 0.008804509063236154, 0.008804509063256116, -0.9691150024934266, -0.5568479179033298, -0.5568479179033498, -0.9691150024933873, 0.4024483358267018, 0.6382186696145582  …  0.5419365898655986, 0.7529730203652125, 0.7529730203651925, 0.5419365898656485, 0.2774134060711553, -0.2774134060711654, -0.47064923446256907, -0.4706492344625691, 0.4706492344625691, 0.4706492344625691], [-0.1408721450121023, -0.14087214501209985, -0.1408721450121023, -0.14087214501209985, -0.13333333333332764, -0.143672232115815, -0.1436722321158125, -0.13333333333333264, -0.2666666666666653, -0.010338898782480455  …  0.09049890913642401, 0.13554687099016938, 0.13554687099016688, 0.09049890913643151, 0.13554687099016688, 0.13554687099016935, 0.3332906247936249, 0.2362763747790916, 0.4303048748081582, 0.3332906247936249])

#@btime Enzyme.gradient(Reverse, $calculate_area, $s, $xpoints, $ypoints, $zpoints)
#   88.625 μs (510 allocations: 59.20 KiB)
