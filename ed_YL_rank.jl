using LinearAlgebra,LinearMaps
using SparseArrays
using Arpack
using ArnoldiMethod
using KrylovKit
using BenchmarkTools
using JLD2

const MyInt = Int64
const MyFloat = Float64

# const L = 8
# const dataPath = "data_rank/"

const L = parse(Int64, ARGS[1])
const dataPath = "/lustre/work/yinghsuan.lin/ed/data_qr/" # NOTE If on cluster set to scratch space

#=
Save/load hard disk to reduce memory usage when measuring ρ.
=#

const dataPathL = dataPath * string(L) * "/"
if !ispath(dataPath)
	mkdir(dataPath)
end
if !ispath(dataPathL)
	mkdir(dataPathL)
end
if !ispath(dataPathL * "prep/")
	mkdir(dataPathL * "prep/")
end

println()
flush(stdout)
println("rank of L=", L)
println()
flush(stdout)

println("available number of threads: ", Threads.nthreads())
println()
flush(stdout)

#=
NOTE To facilitate the code, I always specify the variable types for functions.

The following types are used:
	Bool for anything binary contained in a vector.
		e.g. fusion flag where 0 indicates allowed and 1 disallowed.
	Int8 for anything constant or linear in L contained in a vector.
		e.g. i, pos.
	Int64 for everything else, including variables exponential in L.
		e.g. state, ind, flag, edgeState.
		edgeState assigns one byte for each site, so Int128 needed if L>21.
=#

#=
As Julia's array is 1-based, I use

	"ind" = 1 ... 4^L

as the array index.

"state" is what you obtain by encoding vertex labels using the following mapping,
and takes the values 0 ...  4^L-1.
=#

# y, z, p, m are x, 0, +, - but with ρ draped below at that vertex
const mapping=Dict('x'=>0, '0'=>1, '+'=>2, '-'=>3, 'y'=>0, 'z'=>1, 'p'=>2, 'm'=>3)
const revmapping=Dict(0=>'x', 1=>'0', 2=>'+', 3=>'-')
const startMapping = Dict('0'=>0, '1'=>1, '2'=>2)

#=
The conversion between "ind" and "state" is basically done by considering modulo 4^L.

When L is even, however, some additional care is needed,
since "xxx...x" can come with two "start labels".
My convention is to use
	ind == 4^L  for "x....x" with ρ as the start label
	ind == 1  for "x....x" with 1 as the start label

NOTE For use in zipper, state has two extra bits to encode the start label,
and ~log(L) extra bits to encode the position where ρ is draped below.
	e.g. "x....x" with aρ as the start label and where third x is draped below
		has state with bitstring 0...011010....0 meaning
		0...011   01   0....0
		draped  start  x....x

Now state has two meanings, either with the start/draped info or without.
TODO It is better to use different variable names to distinguish the two.
NOTE I have also decided to change the relation between ind and state to
	ind = state + 1
so that it is less awkward incorporating the extra bits for start/draped.
=#

function stateFromString(s::String,L::Int64=L)
	if length(s)!=L && length(s)!=L+2
		error("length doesn't match")
	end
	if length(s)==L
		s = s*"_0"
	end
	start=startMapping[s[L+2]]
	a=start<<(2*L)
	for i in L : -1 : 1
		if s[i] in ['y','z','p','m']
			a+=(i<<(2*(L+1)))
		end
		a+=(mapping[s[i]]<<(2*(i-1)))
	end
	return a
end

function stringFromState(state::Int64,L::Int64=L)
	below = (state>>(2*(L+1)))
	start = (state>>(2*L)) & 3
	if iseven(L) && (state&(4^L-1))==1
		state -= 1
	end
	s=""
	for i in 1 : L
		if i == below
			if (state&3) == 0
				s*='y'
			elseif (state&3) == 1
				s*='z'
			elseif (state&3) == 2
				s*='p'
			elseif (state&3) == 3
				s*='m'
			else
				error("no match")
			end
		else
			s*=revmapping[(state&3)]
		end
		state>>=2
	end
	s*='_'
	s*=string(start)
	return s
end


#=
Divide-and-conquer to construct basis.
=#

# key is (L, evenxs_left, evenxs_right, Z3 charge)
# value is state (not ind)
basisLego_ = Dict{Tuple{Int8,Bool,Bool,Int8},Vector{Int64}}()

# initialize divide and conquer
for el = false : true
	for er = false : true
		for q = 0 : 2
			basisLego_[(0, el, er, q)] = []
			basisLego_[(1, el, er, q)] = []
		end
	end
end
basisLego_[(1, true, true, 0)] = [ 1 ]
basisLego_[(1, true, true, 1)] = [ 2 ]
basisLego_[(1, true, true, 2)] = [ 3 ]

function setBasisLego!(basisLego::Dict{Tuple{Int8,Bool,Bool,Int8},Vector{Int64}},
	L::Int64, evenxs_left::Bool, evenxs_right::Bool, q::Int64)
	basisLego[(L, evenxs_left, evenxs_right, q)] = []
	if evenxs_right
		# append X
		append!( basisLego[(L, evenxs_left, evenxs_right, q)], getBasisLego(basisLego, L-1, evenxs_left, false, q) )
		# append non-X
		append!( basisLego[(L, evenxs_left, evenxs_right, q)], [ x + (1 << (2*(L-1))) for x in getBasisLego(basisLego, L-1, evenxs_left, true, q) ] )
		append!( basisLego[(L, evenxs_left, evenxs_right, q)], [ x + (2 << (2*(L-1))) for x in getBasisLego(basisLego, L-1, evenxs_left, true, (q+2)%3) ] )
		append!( basisLego[(L, evenxs_left, evenxs_right, q)], [ x + (3 << (2*(L-1))) for x in getBasisLego(basisLego, L-1, evenxs_left, true, (q+1)%3) ] )
	else
		append!( basisLego[(L, evenxs_left, evenxs_right, q)], getBasisLego(basisLego, L-1, evenxs_left, true, q) )
	end
	if ((iseven(L) && !evenxs_left) || (isodd(L) && evenxs_left)) && evenxs_right
		# append non-X to XX...X
		append!( basisLego[(L, evenxs_left, evenxs_right, q)], [ (q+1) << (2*(L-1)) ] )
	end
end

# basis of states (not inds)
function getBasisLego(basisLego::Dict{Tuple{Int8,Bool,Bool,Int8},Vector{Int64}},
	L::Int64, evenxs_left::Bool, evenxs_right::Bool, q::Int64)::Vector{Int64}
	if !haskey(basisLego, (L, evenxs_left, evenxs_right, q))
		setBasisLego!(basisLego, L, evenxs_left, evenxs_right, q)
	end
	return basisLego[(L, evenxs_left, evenxs_right, q)]
end

# basis of states (not inds)
function getBasisLego(basisLego::Dict{Tuple{Int8,Bool,Bool,Int8},Vector{Int64}}, L::Int64)::Vector{Int64}
	res = []
	append!(res, getBasisLego(basisLego, L, true, true, 0))
	append!(res, getBasisLego(basisLego, L, false, false, 0))
	if iseven(L)
		append!(res, [0, 1])
	end
	sort!(res)
	return res
end

const basisPath = dataPathL * "basis.jld2"
if ispath(basisPath)
	println("load basis...")
	flush(stdout)
	@time @load basisPath basis len fromInd
else
	println("compute basis...")
	flush(stdout)
	@time const basis = [ x+1 for x in getBasisLego(basisLego_, L) ]
	const len = length(basis)
	@time const fromInd = Dict((basis[x],x) for x in 1:len)
	@time @save basisPath basis len fromInd
end
println()
flush(stdout)


#=
Construct Hamiltonian.
=#

function trailingXs(state::Int64,L::Int64=L)
	state = state & (4^L-1)
	i=0
	while state!=0
		state>>=2
		i+=1
	end
	return L-i
end

flag_ = zeros(MyInt,len)

function setFlag!(flag::Vector{MyInt},preind::Int64,L::Int64=L)
	ind = basis[preind]
	state=ind-1
	below=(state>>(2*(L+1)))
	start=(state>>(2*L))&3
	state=state&(4^L-1)
	if (start==3) || (below>=L)
		return
	end
	if state==0 && isodd(L)
		return
	end
	evenxs=iseven(trailingXs(state,L))
	tot=start
	if(state==1 && iseven(L))
		state=0
		evenxs=false
	end
	for pos = 0 : L-1
		tot%=3
		a=(state >> (2*pos)) & 3
		if a==0
			if(evenxs)
				flag[preind] |= 1<<pos
			end
			evenxs = ! evenxs
			if ((pos+1)==below || (below!=0 && pos==L-1))
				tot = 3-tot
			end
		else
			if(!evenxs)
				return
			end
			if a==2
				tot+=1
			elseif a==3
				tot+=2
			end
		end
	end
	tot=tot-start
	if tot<0
		tot+=3
	end
	tot%=3
	if state==0 && isodd(L)
		return
	end
end

diag_ = zeros(MyFloat,len)

const ζ = (√(13)+3)/2
const ξ = 1/√ζ
const x = (2-√(13))/3
const z = (1+√(13))/6
const y1 = (5-√(13) - √(6+6*√(13)))/12
const y2 = (5-√(13) + √(6+6*√(13)))/12

const sXX=(0,0)
const sX0=(0,1)
const sXP=(0,2)
const sXM=(0,3)
const s0X=(1,0)
const s00=(1,1)
const s0P=(1,2)
const s0M=(1,3)
const sPX=(2,0)
const sP0=(2,1)
const sPP=(2,2)
const sPM=(2,3)
const sMX=(3,0)
const sM0=(3,1)
const sMP=(3,2)
const sMM=(3,3)

function stateFromInd(ind::Int64,L::Int64=L)
	state=(ind-1)&(4^L-1)
	if state==1 && iseven(L)
		state=0
	end
	return state
end

function nextSite(i::Int64,L::Int64=L)
	j=i+1
	if j==L+1
		j=1
	end
	return j
end

function localStatePair(state::Int64,i::Int64,L::Int64=L)
	j = nextSite(i,L)
	a = (state >> (2*(i-1))) & 3
	b = (state >> (2*(j-1))) & 3
	return a,b
end

function isρ1ρ(flag::Vector{MyInt},preind::Int64,i::Int64)
	return ((flag[preind] >> (i-1)) & 1) == 1
end

function computeDiag!(diag::Vector{MyFloat},flag::Vector{MyInt},preind::Int64)
	ind = basis[preind]
	state=stateFromInd(ind)
	diag[preind]=0
	for i = 1 : L
		sp=localStatePair(state,i)
		if sp==sX0
			diag[preind] -= 1
		elseif sp==s0X
			diag[preind] -= 1
		elseif sp==sXX && isρ1ρ(flag,preind,i)
			diag[preind] -= 1/ζ
		elseif sp==sPM
			diag[preind] -= y1 * y1
		elseif sp==sMP
			diag[preind] -= y2 * y2
		elseif sp==s00
			diag[preind] -= x * x
		elseif sp==s0P
			diag[preind] -= y1 * y1
		elseif sp==sP0
			diag[preind] -= y2 * y2
		elseif sp==sMM
			diag[preind] -= z * z
		elseif sp==s0M
			diag[preind] -= y2 * y2
		elseif sp==sM0
			diag[preind] -= y1 * y1
		elseif sp==sPP
			diag[preind] -= z * z
		end
	end
end

function newInd(state::Int64,i::Int64,sp::Tuple{Int64, Int64})
	(a,b)=sp
	state &= ~(3<<(2*(i-1)))
	state |= (a<<(2*(i-1)))
	j=nextSite(i)
	state &= ~(3<<(2*(j-1)))
	state |= (b<<(2*(j-1)))
	if(state!=0)
		return 1+state
	end
	if(isodd(i))
		return 1
	else
		return 2
	end
end

const prepPath = dataPathL * "prep/prep.jld2"
if ispath(prepPath)
	println("load flag and diag...")
	flush(stdout)
	@time @load prepPath flag_ diag_
else
	println("compute flag and diag...")
	flush(stdout)
	Threads.@threads for preind = 1 : len
		setFlag!(flag_,preind)
		computeDiag!(diag_,flag_,preind)
	end
	@time @save prepPath flag_ diag_
end
println()
flush(stdout)

newPreind(state,i,sp) = fromInd[newInd(state,i,sp)]

TT = Union{Vector{MyFloat},Vector{Float16}}

function sortAndAppendColumn!(
	col::Vector{MyInt},
	row::Vector{MyInt},
	val::TT,
	miniRow::Vector{Int64},
	miniVal::TT
	)
	perm = sortperm(miniRow)
	miniRow = miniRow[perm]
	miniVal = miniVal[perm]
	newMiniRow = MyInt[]
	newMinival = MyFloat[]
	oldr = 0
	v = 0
	cnt = 0
	for r in miniRow
		cnt += 1
		if r == oldr || oldr == 0
			v += miniVal[cnt]
		else
			push!(newMiniRow, oldr)
			push!(newMinival, v)
			v = miniVal[cnt]
		end
		oldr = r
	end
	push!(newMiniRow, oldr)
	push!(newMinival, v)
	push!(col, size(newMiniRow, 1))
	append!(row, newMiniRow)
	append!(val, newMinival)
end

function buildH(diag,flag)
	res = sparse(MyInt[],MyInt[],MyFloat[],len,len)
	col=MyInt[]
	row=MyInt[]
	val=MyFloat[]
	ncol = 1
	for preind = 1 : len
		ind = basis[preind]
		state=stateFromInd(ind)
		# cannot directly append into row, val because row needs to be ordered
		# define miniRow/val and append to row/val after sorting
		miniRow = [preind]
		miniVal = [diag[preind]]
		for i = 1 : L
			sp=localStatePair(state,i)
			if sp==sXX  && isρ1ρ(flag,preind,i)
				append!(miniRow,map(s->newPreind(state,i,s),[sPM,sMP,s00]))
				append!(miniVal,-ξ .* [y1,y2,x])
			elseif sp==sPM
				append!(miniRow,map(s->newPreind(state,i,s),[sXX,sMP,s00]))
				append!(miniVal,-y1 .* [ξ,y2,x])
			elseif sp==sMP
				append!(miniRow,map(s->newPreind(state,i,s),[sXX,sPM,s00]))
				append!(miniVal,-y2 .* [ξ,y1,x])
			elseif sp==s00
				append!(miniRow,map(s->newPreind(state,i,s),[sXX,sPM,sMP]))
				append!(miniVal,-x .* [ξ,y1,y2])
			elseif sp==s0P
				append!(miniRow,map(s->newPreind(state,i,s),[sP0,sMM]))
				append!(miniVal,-y1 .* [y2,z])
			elseif sp==sP0
				append!(miniRow,map(s->newPreind(state,i,s),[s0P,sMM]))
				append!(miniVal,-y2 .* [y1,z])
			elseif sp==sMM
				append!(miniRow,map(s->newPreind(state,i,s),[s0P,sP0]))
				append!(miniVal,-z .* [y1,y2])
			elseif sp==s0M
				append!(miniRow,map(s->newPreind(state,i,s),[sM0,sPP]))
				append!(miniVal,-y2 .* [y1,z])
			elseif sp==sM0
				append!(miniRow,map(s->newPreind(state,i,s),[s0M,sPP]))
				append!(miniVal,-y1 .* [y2,z])
			elseif sp==sPP
				append!(miniRow,map(s->newPreind(state,i,s),[s0M,sM0]))
				append!(miniVal,-z .* [y2,y1])
			end
		end

		sortAndAppendColumn!(col, row, val, miniRow, miniVal)

		if (preind % (len / 10)) == 1 || preind == len
			num = res.colptr[ncol]
			for c in col
				ncol += 1
				num += c
				res.colptr[ncol] = num
			end
			append!(res.rowval, row)
			append!(res.nzval, val)
			col=MyInt[]
			row=MyInt[]
			val=MyFloat[]
		end
		append!(res.rowval, row)
		append!(res.nzval, val)
		row=MyInt[]
		val=MyFloat[]
	end
	return res
end

function Hfunc!(C,B,diag::Vector{MyFloat},flag::Vector{MyInt})
	Threads.@threads for preind = 1 : len
		C[preind] = diag[preind] * B[preind]
	end
	for i = 1 : L
		Threads.@threads for preind = 1 : len
			ind = basis[preind]
			state=stateFromInd(ind)
			sp=localStatePair(state,i)
			if sp==sXX  && isρ1ρ(flag,preind,i)
				C[preind]-= ξ * (
				     y1 * B[newPreind(state,i,sPM)] +
					 y2 * B[newPreind(state,i,sMP)] +
					 x  * B[newPreind(state,i,s00)]
			    )
			elseif sp==sPM
				C[preind]-= y1*(
					 ξ * B[newPreind(state,i,sXX)] +
					 y2* B[newPreind(state,i,sMP)] +
					 x * B[newPreind(state,i,s00)]
			    )
			elseif sp==sMP
				C[preind]-= y2*(
					 ξ * B[newPreind(state,i,sXX)] +
					 y1* B[newPreind(state,i,sPM)] +
					 x * B[newPreind(state,i,s00)]
			    )
			elseif sp==s00
				C[preind]-= x*(
					 ξ * B[newPreind(state,i,sXX)] +
					 y1* B[newPreind(state,i,sPM)] +
					 y2* B[newPreind(state,i,sMP)]
			    )
			elseif sp==s0P
				C[preind]-= y1*(
					y2 * B[newPreind(state,i,sP0)] +
					z  * B[newPreind(state,i,sMM)]
				)
			elseif sp==sP0
				C[preind]-= y2*(
					y1 * B[newPreind(state,i,s0P)] +
					z  * B[newPreind(state,i,sMM)]
				)
			elseif sp==sMM
				C[preind]-= z*(
					y1 * B[newPreind(state,i,s0P)] +
					y2 * B[newPreind(state,i,sP0)]
				)
			elseif sp==s0M
				C[preind]-= y2*(
					y1 * B[newPreind(state,i,sM0)] +
					z  * B[newPreind(state,i,sPP)]
				)
			elseif sp==sM0
				C[preind]-= y1*(
					y2 * B[newPreind(state,i,s0M)] +
					z  * B[newPreind(state,i,sPP)]
				)
			elseif sp==sPP
				C[preind]-= z*(
					y2 * B[newPreind(state,i,s0M)] +
					y1 * B[newPreind(state,i,sM0)]
				)
			end
		end
	end
end

HPath = dataPathL * "H.jld2"
if ispath(HPath)
	println("load H...")
	flush(stdout)
	@load HPath H
else
	println("build H...")
	flush(stdout)
	@time H=buildH(diag_,flag_)
	@save HPath H
end
println()
flush(stdout)
H = -H # ferro

flush(stdout)
@time println("Number of ground states from rank: ", size(H,1) - rank(H))
println()
flush(stdout)
@time println("Number of ground states from rank qr: ", size(H,1) - rank(qr(H)))
println()
flush(stdout)
