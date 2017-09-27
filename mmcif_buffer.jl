const nucleic = ["GUN","ADE","CYT","THY","URA","DA","DC","DG","DT","DU","A","C","G","T","U",
                "AMP","ADP","ATP","CDP","CTP","GMP","GDP","GTP","TMP","TTP","UMP","UDP","UTP"]
const water = ["HOH","H2O"]
const aas = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
const lipids = []
#and more

function isheader(linestr::String)
    @inbounds return linestr[1] == '_' || linestr[1] == '#'
end

veclength2(ax,bx,ay,by,az,bz) = (ax-bx)^2 + (ay-by)^2 + (az-bz)^2

const _default_delims = ' '
function splitmmcif(str::AbstractString,atomid,resn,chain,resid,inscode,xlist,ylist,zlist,strs,coords)
    #strs = fill("",5)
    #coords = zeros(Float64,3)
    splitter = _default_delims

    i = start(str)
    n = endof(str)
    r = search(str,splitter,i)
    j, k = first(r), nextind(str,last(r))
    colno = 1

    while 0 < j <= n
        if i < k
            if i < j
                if colno == atomid
                    strs[1] = str[i:j-1]
                elseif colno == resn
                    strs[2] = str[i:j-1]
                elseif colno == chain
                    strs[3] = str[i:j-1]
                elseif colno == resid
                    strs[4] = str[i:j-1]
                elseif colno == inscode
                    strs[5] = str[i:j-1]
                elseif colno == xlist
                    coords[1] = coordsfloat(str[i:j-1])
                elseif colno == ylist
                    coords[2] = coordsfloat(str[i:j-1])
                elseif colno == zlist
                    coords[3] = coordsfloat(str[i:j-1])
                end
                colno += 1
            end
            i = k
        end
        (k <= j) && (k = nextind(str,j))
        r = search(str,splitter,k)
        j, k = first(r), nextind(str,last(r))
    end
    if !done(str,i)
        if colno == atomid
            strs[1] = str[i:j-1]
        elseif colno == resn
            strs[2] = str[i:j-1]
        elseif colno == chain
            strs[3] = str[i:j-1]
        elseif colno == resid
            strs[4] = str[i:j-1]
        elseif colno == inscode
            strs[5] = str[i:j-1]
        elseif colno == xlist
            coords[1] = coordsfloat(str[i:j-1])
        elseif colno == ylist
            coords[2] = coordsfloat(str[i:j-1])
        elseif colno == zlist
            coords[3] = coordsfloat(str[i:j-1])
        end
    end
    #return strs,coords
end

function coordsfloat(str)
    ptr = pointer(str)
    #parse str(xxxxx.xxxxx) -> float
    negaflag = unsafe_load(ptr,1) == 0x2d ? -1 : 1
    n = endof(str)
    i = (3 - negaflag) >> 1
    decspace = 0
    isdot = false
    val = 0.0

    for idx = i:n
        digit = unsafe_load(ptr,idx)
        if digit == 0x2e
            decspace = 0
            isdot = true
            continue
        end
        val *= 10
        val += digit - 0x30
        decspace += 1
    end
    if isdot
        val /= 10^decspace
    end

    val *= negaflag
    return val
end


#--------------------------------------MAIN
function read_mmcif(datalines,thresholdlength=5.0::Float64,targetmoleculelists=nucleic) #change here or specify the target explicitly
    const thresholdlengthpow2 = thresholdlength^2

    #find the _atom_site-loop and memo column no of coords
    const atomcolumn = findfirst(x -> contains(x,"_atom_site."),datalines) -1 #line of "loop_"
    const datastart = findnext(!isheader,datalines,atomcolumn+1)
    const datalast = findnext(isheader,datalines,datastart) -1 #previous line of "#"

    const atomid = findnext(x -> contains(x,"_atom_site.auth_atom_id"),datalines,atomcolumn) - atomcolumn
    const resn = findnext(x -> contains(x,"_atom_site.auth_comp_id"),datalines,atomcolumn) - atomcolumn
    const chain = findnext(x -> contains(x,"_atom_site.auth_asym_id"),datalines,atomcolumn) - atomcolumn
    const resid = findnext(x -> contains(x,"_atom_site.auth_seq_id"),datalines,atomcolumn) - atomcolumn
    const inscode = findnext(x -> contains(x,"_atom_site.pdbx_PDB_ins_code"),datalines,atomcolumn) - atomcolumn
    const xlist = findnext(x -> contains(x,"_atom_site.Cartn_x"),datalines,atomcolumn) - atomcolumn
    const ylist = findnext(x -> contains(x,"_atom_site.Cartn_y"),datalines,atomcolumn) - atomcolumn
    const zlist = findnext(x -> contains(x,"_atom_site.Cartn_z"),datalines,atomcolumn) - atomcolumn

    #prepare array to store (textbase/coordinary) data of (amino acids/ligands)
    aalist = fill("",(5,(datalast-datastart)))
    aacoords = Array{Float64,2}(3,(datalast-datastart))
    liglist = fill("",(5,(datalast-datastart)))
    ligcoords = Array{Float64,2}(3,(datalast-datastart))

    aacount = 1
    ligcount = 1
    aap = pointer(aacoords)
    ligp = pointer(ligcoords)
    strpart = fill("",5)
    coordspart = zeros(Float64,3)

    for lineno = datastart:datalast
        #split one line text to array
        splitmmcif(datalines[lineno],atomid,resn,chain,resid,inscode,xlist,ylist,zlist,strpart,coordspart)

        if strpart[2] in aas #canonical aa.
            #I can't use pointer for aa/lig-list, because something cause error.
            @inbounds aalist[:,aacount] = strpart
            @inbounds aacoords[:,aacount] = coordspart
            aacount += 1
        end
        if strpart[2] in targetmoleculelists
            @inbounds liglist[:,ligcount] = strpart
            @inbounds ligcoords[:,ligcount] = coordspart
            ligcount += 1
        end
    end

    aacount -= 1
    ligcount -= 1

    #cut unnecessary region of coords matrix.
    #text-base information matrixes may not be needed to cut, so I don't.
    @inbounds aacoords = aacoords[:,1:aacount]
    @inbounds ligcoords = ligcoords[:,1:ligcount]


    #-------------------------------------------O(N) search from here
    @inbounds xmin = min(minimum(aacoords[1,:]),minimum(ligcoords[1,:]))
    @inbounds xmax = max(maximum(aacoords[1,:]),maximum(ligcoords[1,:]))
    @inbounds ymin = min(minimum(aacoords[2,:]),minimum(ligcoords[2,:]))
    @inbounds ymax = max(maximum(aacoords[2,:]),maximum(ligcoords[2,:]))
    @inbounds zmin = min(minimum(aacoords[3,:]),minimum(ligcoords[3,:]))
    @inbounds zmax = max(maximum(aacoords[3,:]),maximum(ligcoords[3,:]))

    #make the grid, it's like a hash function of 3d-space of the mmcif
    xgrid = ceil(Int,(xmax-xmin)/thresholdlength)
    ygrid = ceil(Int,(ymax-ymin)/thresholdlength)
    zgrid = ceil(Int,(zmax-zmin)/thresholdlength)
    numgrid = xgrid * ygrid * zgrid

    #using sortedbuffer (, count) and  indexes to find of aacoords
    count = zeros(Int64,numgrid) #the no of aa. atoms in each grid
    indexes = zeros(Int64,numgrid+1) #the 1st index of each grid in sortedbuffer
    particlepos = zeros(Int64,aacount) #the grid of each aa. atom
    sortedbuffer = zeros(Int64,aacount) #the list of aa. atoms sorted by grid no.

    for i = 1:aacount
        @inbounds ix = floor(Int, (aacoords[1,i] - xmin) / thresholdlength)
        @inbounds iy = floor(Int, (aacoords[2,i] - ymin) / thresholdlength)
        @inbounds iz = floor(Int, (aacoords[3,i] - zmin) / thresholdlength)

        #grid no. index  "+1": julia is 1-origin.
        idx = floor(Int, ix + iy * xgrid + iz * xgrid * ygrid) + 1
        @inbounds count[idx] += 1
        @inbounds particlepos[i] = idx
    end

    sum = 0
    for i in 1:numgrid
        @inbounds sum += count[i]
        @inbounds indexes[i+1] = sum
    end

    register = zeros(Int64,numgrid) #temp. values how many aa. atoms are already put into sortedbuffer

    #SORT
    for i in 1:aacount
        @inbounds pos = particlepos[i]
        @inbounds j = indexes[pos] + register[pos] + 1
        @inbounds sortedbuffer[j] = i
        @inbounds register[pos] += 1
    end

    #redefine pointers
    aap = pointer(aacoords)
    ligp = pointer(ligcoords)

    #--------------------------------------------------------
    for tgtidx = 1:ligcount

        @inbounds target = liglist[:,tgtidx]
        #@inbounds tgtcoords = ligcoords[:,tgtidx]
        tx = unsafe_load(ligp,3*tgtidx-2)
        ty = unsafe_load(ligp,3*tgtidx-1)
        tz = unsafe_load(ligp,3*tgtidx)

        #find the grid to where the target atom should be belong
        ix = floor(Int, (tx - xmin) / thresholdlength)
        iy = floor(Int, (ty - ymin) / thresholdlength)
        iz = floor(Int, (tz - zmin) / thresholdlength)

        #Because grid length is equal to thresholdlength,
        #the sought aa. atoms within the threshold length must be in the same or the 26-adjacent grid(s).
        for bufx = -1:1
            for bufy = -1:1
                for bufz = -1:1
                    meshsearch(thresholdlengthpow2,aalist,aap,target,tx,ty,tz,sortedbuffer,count,indexes,xgrid,ygrid,zgrid,ix+bufx,iy+bufy,iz+bufz)
                end
            end
        end
    end

end

function meshsearch(thresholdlengthpow2,aalist,aap,target,tx,ty,tz,sortedbuffer,count,indexes,xgrid,ygrid,zgrid,ix,iy,iz)
    if ix < 0 || xgrid <= ix || iy < 0 || ygrid <= iy || iz < 0 || zgrid <= iz
        return 1
    end

    #index of searching grid
    id2 = ix + iy * xgrid + iz * xgrid * ygrid + 1
    @inbounds for m in indexes[id2]+1:indexes[id2]+count[id2]
        @inbounds aaidx = sortedbuffer[m]

        ax = unsafe_load(aap,3*aaidx-2)
        ay = unsafe_load(aap,3*aaidx-1)
        az = unsafe_load(aap,3*aaidx)

        #slight faster than bother to find root
        #lenpow2 = veclength2(vectors)
        lenpow2 = veclength2(ax,tx,ay,ty,az,tz)
        if lenpow2 <= thresholdlengthpow2
            @inbounds aainfo1 = aalist[1,aaidx]
            @inbounds aainfo2 = aalist[2,aaidx]
            @inbounds aainfo3 = aalist[3,aaidx]
            @inbounds aainfo4 = aalist[4,aaidx]
            @inbounds aainfo5 = aalist[5,aaidx]
            #only print 5 decimal space: you can change .[5]f
            @printf "%s.%s.%s.%s.%s:%s.%s.%s.%s:%.5f \n" aainfo2 aainfo3 aainfo5 aainfo4 aainfo1 target[2] target[3] target[4] target[1] (lenpow2^0.5)
        end
    end
end



if length(ARGS) == 1
    readcif = readlines(ARGS[1])
    read_mmcif(readcif)
elseif length(ARGS) == 2
    readcif = readlines(ARGS[1])
    threlen = parse(Float64, ARGS[2])
    read_mmcif(readcif,threlen)
else
    print("Input MMcif file needed: \n /your/path/to/julia mmcif.jl *.cif (threshold_length[Ang.] : optional, default = 5.0)")
end
#=
hoge = readlines("1vy4.cif")
read_mmcif(hoge)
for i=1:5
    @time read_mmcif(hoge)
end
#Profile.print()
=#
