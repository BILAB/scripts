#Julia v0.6.2
#Julia で Pkg.add("PyCall"), python3 でMDAnalysis MDAnalysisTests numpy が必要

using PyCall
@pyimport MDAnalysis as mda
@pyimport MDAnalysis.topology.guessers as mda_tg

#考慮対象外の水とイオン
#他に考慮しない残基名等あればここに追記してください
#最初に回転したりなんかしたりするなら不要？
const wation = ["HOH","H2O","WAT","Na+","Cl-"]

const boltzmann =  1.380649e-23 / 4.184
const boltzmannpermol = 1.380649 * 6.02214076 / 4.184

function dihedral_divide(p1, p2, p3, p4)
    v1 = p1-p2
    v2 = p3-p2
    v3 = p4-p3

    normalize!(v2)
    # v1/3 minus component that aligns with v2
    v1 -= (v1⋅v2)v2
    v3 -= (v3⋅v2)v2

    x = v1⋅v3
    y = (v2×v1)⋅v3

    return atan2(y,x) #(-π,π]
end
dihedral(p::Array{Float64,2}) = dihedral_divide(p[1,:], p[2,:], p[3,:], p[4,:])
dihedral(p::Array{Float64,2}, i1::Int, i2::Int, i3::Int, i4::Int) = dihedral_divide(p[i1,:], p[i2,:], p[i3,:], p[i4,:])
dihedral(p::Array{Float64,2}, iv::Vector{Int}) = dihedral(p, iv[1], iv[2], iv[3], iv[4])

function angle_divide(p1, p2, p3)
    v1 = p2-p1
    v2 = p2-p3
    normalize!(v1)
    normalize!(v2)
    return acos(v1⋅v2) #[0,π]
end
angle(p::Array{Float64,2}) = angle_divide(p[1,:], p[2,:], p[3,:])
angle(p::Array{Float64,2}, i1::Int, i2::Int, i3::Int) = angle_divide(p[i1,:], p[i2,:], p[i3,:])
angle(p::Array{Float64,2}, iv::Vector{Int}) = angle(p, iv[1], iv[2], iv[3])

function improper_divide(p1, p2, p3, p4)
    # p1 = root, p2 = hub, p3 & p4 = leaf
    # dihedral of plane(p1,p2,p3) and plane(p1,p2,p4)

    v1 = p2-p1
    v2 = p3-p2
    v3 = p4-p2

    normalize!(v1)
    v2 -= (v1⋅v2)v1
    v3 -= (v1⋅v3)v1

    x = v2⋅v3
    y = (v1×v2)⋅v3

    return atan2(y,x)
end
improper(p::Array{Float64,2}) = improper_divide(p[1,:], p[2,:], p[3,:], p[4,:])
improper(p::Array{Float64,2}, i1::Int, i2::Int, i3::Int, i4::Int) = improper_divide(p[i1,:], p[i2,:], p[i3,:], p[i4,:])
improper(p::Array{Float64,2}, iv::Vector{Int}) = improper(p, iv[1], iv[2], iv[3], iv[4])

bondlen(p::Array{Float64,2}, i1::Int, i2::Int) = vecnorm(p[i1,:] - p[i2,:])
bondlen(p::Array{Float64,2}, iv::Vector{Int}) = bondlen(p, iv[1], iv[2])

function connected_hydrogen(namelist::Vector{String}, connect, parent)
    ret = Int[]
    for i in connect[parent]
        if namelist[i][1] == 'H'
            push!(ret,i)
        end
    end
    return ret
end

const nm_conv = log(10.0^-9)
function jacobian(tm_val, angls_idx)
    tmpjac = 0.0
    #tmpjac += 2 * (log(tm_val[1]) )
    for i in 1:angls_idx
        #tmpjac += 2 * (log(tm_val[i+1])) 
        tmpjac += log( sin(tm_val[i+angls_idx+1]) )
    end
    return tmpjac
end

function load_coord!(vals, io::IO, n_trajectory::Int, head::Int)
    seekstart(io)
    skip(io, 8 * (head-1) * n_trajectory )
    read!(io, vals)
end

function shannon(vals::Vector{Float64}, bins::Int, pitch::T where T<:AbstractFloat)
    nval = length(vals)
    ret = 0.0
    ctr = zeros(Int,bins+1)
    for i in vals
        ctr[floor(Int,i/pitch) + 1] += 1
    end
    for i in ctr
        if i == 0; continue; end
        p_i = Float64(i) / nval
        ret += p_i * log(p_i)
    end
    return -ret
end

function shannon_bond(vals::Vector{Float64}, bins::Int)
    vals .-= minimum(vals)
    pitch = maximum(vals) / bins
    return shannon(vals, bins, pitch)
end

function shannon_angle(vals::Vector{Float64}, bins::Int)
    vals .-= minimum(vals)
    pitch = maximum(vals) / bins
    return shannon(vals, bins, pitch)
end

function shannon_dihedral(vals::Vector{Float64}, bins::Int)
    vals .+= π
    pitch = 2π / bins
    return shannon(vals, bins, pitch)
end

function fast_cor(n_trajectory::Int, idx::Int, idy::Int)

    sx = unsafe_load(pmean,idx)
    sy = unsafe_load(pmean,idy)

    xy = 0.0

    @simd for i = 1:n_trajectory
        xi = unsafe_load(px,i)
        yi = unsafe_load(py,i)
        xi -= sx
        yi -= sy
        #xx += xi * xi
        #yy += yi * yi
        xy += xi * yi
    end

    xx = unsafe_load(psp,idx)
    yy = unsafe_load(psp,idy)

    return (xy / (sqrt(xx*yy)))
end

ff14types = Dict(
       "CA" => 1.5,
       "N3" => 1.1,
       "CN" => 1.5,
       "CT" => 1.5,
       "CO" => 1.5,
       "CB" => 1.5,
       "HA" => 0.4,
       "HC" => 0.4,
       "H4" => 0.4,
       "HO" => 0.4,
       "OH" => 1.05,
       "CX" => 1.5,
       "O2" => 1.05,
       "NA" => 1.1,
       "H1" => 0.4,
       "HP" => 0.4,
       "3C" => 1.5,
       "C8" => 1.5,
       "C*" => 1.5,
       "2C" => 1.5,
       "CW" => 1.5,
       "N2" => 1.1)

#テストのお約束
#PSF = mdad.PSF #data/adk.psf
#DCD = mdad.DCD #data/adk_dims.dcd
TPR = ARGS[1] 
TRR = ARGS[2] 
const test = mda.Universe(TPR,TRR)

const atomname = convert.(String, test[:atoms][:names]) #原子名
const resname = convert.(String, test[:atoms][:resnames]) #原子ごとの残基名
const n_atoms = length(atomname)

const allbonds = ARGS[4] == "0" ? test[:bonds] : mda_tg.guess_bonds(test[:atoms],test[:atoms][:positions],vdwradii=ff14types)

const n_bonds = length(allbonds)
@show const n_trajectory = length(test[:trajectory])
const temperature = float(ARGS[3]) 



#結合情報の格納　グラフの次数が高々4なのでリストで持つ
connect = [Int[] for i=1:n_atoms]
@inbounds for bond in allbonds
    bondno = bond[:atoms][:ix]
    push!(connect[bondno[1] + 1], bondno[2] + 1) #python : 0-index julia : 1-index
    push!(connect[bondno[2] + 1], bondno[1] + 1)
end

#参照先の格納
bonds = zeros(Int,(n_atoms-1,2))
angls = zeros(Int,(n_atoms-2,3))
dihds = zeros(Int,(n_atoms-3,4))
bonds_idx = angls_idx = dihds_idx = 1

#木のrootについて ad hoc な取扱い。
#グラフの次数で見て、0ならばイオンか何かなので飛ばす
#1ならばそれをrootに、4ならば多分N末N原子なので隣にある水素をどれか取って無理やりrootに
#タンパク質が存在しないと困るのでグラフの次数2-3はエラーを吐く　←改善の余地あり
rootidx = 1
seeking = true
while seeking
    if length(connect[rootidx]) == 4
        for i in connect[rootidx]
            if length(connect[i]) == 1
                rootidx = i
                seeking = false
                break
            end
        end
    elseif length(connect[rootidx]) == 1
        seeking = false
    elseif 2 <= length(connect[rootidx]) <= 3
        println( STDERR, "invalid connection order $(connect[rootidx]) in initializing the internal coordinate tree")
        return 1
    end
    rootidx += 1
end
rootidx -= 1

#rootからBFSして内部座標の構築
#rootを次数0として、次数1のノードは1つになるように取った
#グラフの次数2のノードが複数あるとdihedralのとり方が非自明になるので、どこか一個固定してimproperとしてみる
two_top = 0
#親
parent = zeros(Int,n_atoms)
parent[rootidx] = -1
#juliaだと特別dequeとかあるわけではなくvectorで全部ごり押すらしい
queue = [rootidx]
#念の為スコープの外で宣言
top = rootidx

finished = false
@inbounds while !finished
    if length(queue) != 0
        while length(queue) > 0
            top = shift!(queue)

            if parent[top] != -1
                bonds[bonds_idx,1] = top
                bonds[bonds_idx,2] = parent[top]
                bonds_idx += 1
            end
            if bonds_idx == 3
                two_top = top
            end

            if bonds_idx >= 3
                angls[angls_idx,1] = top
                angls[angls_idx,2] = parent[top]
                angls[angls_idx,3] = parent[parent[top]]
                angls_idx += 1
            end

            if bonds_idx >= 4

                #メチル基判定で使う
                #親ノードの接続している原子のうち原子名の1文字目が'H'であるインデックスのリスト
                #note:空になりうるので長さ判定を先に入れること
                cnh = connected_hydrogen(atomname, connect, parent[top])
                #次数2判定
                if parent[parent[parent[top]]] == -1
                    dihds[dihds_idx,1] = parent[parent[top]]
                    dihds[dihds_idx,2] = parent[top]
                    dihds[dihds_idx,3] = two_top
                    dihds[dihds_idx,4] = top
                #メチル基Hの判定　improperになるのはメチル基Hの2番手以降
                elseif atomname[top][1] == 'H' && atomname[parent[top]][1] == 'C' && length(cnh) == 3 && cnh[1] != top
                    dihds[dihds_idx,1] = parent[parent[top]]
                    dihds[dihds_idx,2] = parent[top]
                    dihds[dihds_idx,3] = cnh[1]
                    dihds[dihds_idx,4] = top
                #それ以外は通常のtorsion
                else
                    dihds[dihds_idx,1] = parent[parent[parent[top]]]
                    dihds[dihds_idx,2] = parent[parent[top]]
                    dihds[dihds_idx,3] = parent[top]
                    dihds[dihds_idx,4] = top
                end

                dihds_idx += 1
            end #add dihedral

            for i in connect[top]
                if parent[i] == 0
                    parent[i] = top
                    push!(queue,i)
                end
            end

        end #BFS while
    else # = queue is empty

        #チェインの切れ目とかが存在した場合、直前のところに無理やりつなぐ
        for i in 1:n_atoms
            if parent[i] == 0 && !(resname[i] in wation)
                parent[i] = top
                finished = !finished
                break
            end
        end
        #最後まで探索してしまった場合1回しかflipされないのでそのまま最外ループを抜ける
        finished = !finished

    end # if queue is not empty
end #while !finished

bonds_idx -= 1
angls_idx -= 1
dihds_idx -= 1
const n_coords = bonds_idx + angls_idx + dihds_idx
@show const n_nowation = (n_coords+6)/3

#S_{ind} の項
#各座標の時間発展を全部取ると確実にメモリからはみ出るので
#一度中間ファイルに全部書いてから一般化座標の一個ずつ足す
#水はないから多少はファイルサイズ減るはず
#-kBは括りだして後で掛ける

#S_{jac} の項 : jacobianがかかってくる項
#b2^2 ∏_{i=3}^{N} b_i^2 sin θ_i がjacobianで、それにlnがかかっているのでlnとって総和
#あとでpdf考慮するためにn_trajectoryで割ってkB掛ける

S_ind = 0.0
S_jac = 0.0

tempfile = open("mdentropy_tmp","w+") #w,r,create,truncate
unsafe_write(tempfile,Ref(UInt8(0)), 8 * n_coords * n_trajectory) #確保

tm_val = zeros(Float64, n_coords) #一般化座標の格納先
@inbounds for tm = 1:n_trajectory
    dump = test[:trajectory][tm] #フレーム移動のお約束

    coords = Float64.(test[:atoms][:positions])

    for i in 1:bonds_idx
        tm_val[i] = bondlen(coords,bonds[i,:])
    end
    for i in 1:angls_idx
        tm_val[i + bonds_idx] = angle(coords,angls[i,:])
    end
    for i in 1:dihds_idx
        i_dihds = dihds[i,:]
        #二面角がproper ⇔ 二面角計算の3,4番が結合している は逆が一般には偽だが
        #BFSで構築した都合上3員環が存在しない仮定を置くと成立する
        if parent[i_dihds[4]] == i_dihds[3]
            tm_val[i + bonds_idx + angls_idx] = dihedral(coords,dihds[i,:])
        else
            tm_val[i + bonds_idx + angls_idx] = improper(coords,dihds[i,:])
        end
    end

    #coord:元 tm_val:内部


    #情報落ち防止のためJ(b,θ)は一旦計算してからS_jacに足す
    S_jac += jacobian(tm_val,angls_idx)

    #各座標ごとにバラして格納
    #magic number の 8 は Float64 = (cのfloat) のbyte数
    seekstart(tempfile)
    skip(tempfile, 8*(tm-1) )
    p = pointer(tm_val)
    for i in 1:n_coords
        unsafe_write(tempfile,p,8)
        skip(tempfile,8*(n_trajectory-1))
        p += 8
    end
end

S_jac /= n_trajectory
S_jac *= boltzmannpermol
println()
println(STDOUT,"S_jac :", S_jac)
#seekstart(tempfile)

#15 = bin数
#あとで変数に変える
vals = zeros(Float64,n_trajectory)
for i = 1:bonds_idx
    load_coord!(vals, tempfile, n_trajectory, i)
    S_ind += shannon_bond(vals,15)
end
for i = 1:angls_idx
    load_coord!(vals, tempfile, n_trajectory, i + bonds_idx)
    S_ind += shannon_angle(vals,15)
end
for i = 1:dihds_idx
    load_coord!(vals, tempfile, n_trajectory, i + bonds_idx + angls_idx)
    S_ind += shannon_dihedral(vals,15)
end
S_ind *= boltzmannpermol
println()
println(STDOUT,"S_ind :", S_ind)


#ここから重いので対策
#繰り返し利用する平均と残差平方和を別に求めておき、相関行列を求めるときに都度持ってくる。

corr_matrix = zeros(Float64,(n_coords,n_coords))
pcorr = pointer(corr_matrix)

#作業スペース
x = zeros(Float64,n_trajectory)
y = zeros(Float64,n_trajectory)
const px = pointer(x)
const py = pointer(y)

const coef = 1.0

#平均
means = zeros(Float64,n_coords)
const pmean = pointer(means)
#平方和
sumpow = zeros(Float64,n_coords)
const psp = pointer(sumpow)

#xを使って上2つの作成
seekstart(tempfile)
for i=1:n_coords
    mx = 0.0
    read!(tempfile, x)
    @simd for j = 1:n_trajectory
        mx += unsafe_load(px,j)
    end
    mx /= n_trajectory
    unsafe_store!(pmean,mx,i)

    psx = 0.0
    @simd for j = 1:n_trajectory
        psx += (unsafe_load(px,j) - mx)^2
    end
    unsafe_store!(psp,psx,i)
end

S_indQH = 0.0
for i=1:n_coords
    @inbounds S_indQH += log(2*pi*e*sumpow[i])
end
S_indQH *= boltzmannpermol/2
println("S_indQH: ",S_indQH)
println()

for i = 1:n_coords
    #corr_matrix[i,i] = coef
    unsafe_store!(pcorr, coef, (i-1)*n_coords+i)
    

    load_coord!(x, tempfile, n_trajectory, i)

    for j=i+1:n_coords
        #一般化座標の連番に並んでいるので、yの大きさだけ読めばj番目の座標データが入る
        read!(tempfile, y)

        #ポインタを明示的に使って高速化
        #tmp = coef * cor(x,y)
        tmp = coef * fast_cor(n_trajectory, i, j)

        #debug
        #@assert abs( tmp - cor(x,y) * coef) < 1e-9 * coef

        #corr_matrix[j,i] = tmp　
        unsafe_store!(pcorr, tmp, (i-1)*n_coords+j)
    end
end

corr_matrix = Symmetric(corr_matrix, :L)
#対称行列なのでBunch-Kauchman分解ができる
logdetcoefcorr, hoge = logabsdet(bkfact(corr_matrix))

logdetcorr = logdetcoefcorr - log(coef)*n_coords

S_corr = boltzmannpermol * logdetcorr / 2
println("S_corr: ",S_corr)

S_P = ((3n_nowation - 6.0) / 2 * log(temperature) + 1.5 * n_nowation ) * boltzmannpermol

println("S_P: ",S_P)

println("---------------------------------")
println("S_Conf: ", S_P+S_ind+S_corr+S_jac)

