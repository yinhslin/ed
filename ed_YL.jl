using LinearAlgebra,LinearMaps
using SparseArrays
using ArnoldiMethod
using Arpack
# using Profile
# using Traceur

const L=9

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
	if length(s)!=L || length(s)!=L+2
		error("length doesn't match")
	end
	if length(s)==L
		s = s*"_0"
	end
	start=Int64(s[L+2])
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
TODO trailingXs is always used to infer whether start label is type 1 or ρ.
Might as well define a function that also checks whether state is 1 when L even.
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

# Changed to Int64 since Int32 will not be enough even for L=15.
flag_ = zeros(Int64,4^L)

#=

Modified to allow for nontrivial start/draped used in the zipper.

But zipper only cares about the main flag, which is computed by setFusionFlag!,
so might as well keep Yuji's original definition.

=#
function setFlag!(flag::Vector{Int64},ind::Int64,L::Int64=L)
	flagshift=L+2
	state=ind-1
	below=(state>>(2*(L+1)))
	start=(state>>(2*L))&3
	state=state&(4^L-1)
	if (start==3) || (below>=L)
		flag[ind] |= 3 << flagshift
		return
	end
	if state==0 && isodd(L)
		flag[ind] |= 3 << flagshift
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
				flag[ind] |= 1<<pos
			end
			evenxs = ! evenxs
			if ((pos+1)==below || (below!=0 && pos==L-1))
				tot = 3-tot
			end
		else
			if(!evenxs)
				# not allowed
				flag[ind] |= 3 << flagshift
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
		flag[ind] |= 3 << flagshift
		return
	end
	flag[ind] |= tot << flagshift
end

diag_ = zeros(Float64,4^L)

const U = 10	# suppression factor for forbidden state
const Z = 10	# suppression factor for states charged under twisted Z3
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

function mainFlag(flag::Vector{Int64},ind::Int64,L::Int64=L)::Int8
	flagshift=L+2
	return (flag[ind] >> flagshift) & 3
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

function isρ1ρ(flag::Vector{Int64},ind::Int64,i::Int64)
	return ((flag[ind] >> (i-1)) & 1) == 1
end

function computeDiag!(diag::Vector{Float64},flag::Vector{Int64},ind::Int64)
	state=stateFromInd(ind)
	fl = mainFlag(flag,ind)
	if fl==3
		diag[ind]=U
		return
	elseif fl==1 || fl==2
		diag[ind]=Z
		return
	end
	diag[ind]=0
	for i = 1 : L
		sp=localStatePair(state,i)
		if sp==sX0
			diag[ind] -= 1
		elseif sp==s0X
			diag[ind] -= 1
		elseif sp==sXX && isρ1ρ(flag,ind,i)
			diag[ind] -= 1/ζ
		elseif sp==sPM
			diag[ind] -= y1 * y1
		elseif sp==sMP
			diag[ind] -= y2 * y2
		elseif sp==s00
			diag[ind] -= x * x
		elseif sp==s0P
			diag[ind] -= y1 * y1
		elseif sp==sP0
			diag[ind] -= y2 * y2
		elseif sp==sMM
			diag[ind] -= z * z
		elseif sp==s0M
			diag[ind] -= y2 * y2
		elseif sp==sM0
			diag[ind] -= y1 * y1
		elseif sp==sPP
			diag[ind] -= z * z
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

const Hairetsu=SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}

function Hfunc!(C::Hairetsu,B::Hairetsu,diag::Vector{Float64},flag::Vector{Int64})
	for ind = 1 : 4^L
		C[ind] = diag[ind] * B[ind]
	end
	for ind = 1 : 4^L
		if  mainFlag(flag,ind) !=0
			continue
		end
		state=stateFromInd(ind)
		for i = 1 : L
			sp=localStatePair(state,i)
			if sp==sXX  && isρ1ρ(flag,ind,i)
				C[newInd(state,i,sPM)] -= ξ * y1 * B[ind]
				C[newInd(state,i,sMP)] -= ξ * y2 * B[ind]
				C[newInd(state,i,s00)] -= ξ * x * B[ind]
			elseif sp==sPM
				C[newInd(state,i,sXX)] -= y1 * ξ * B[ind]
				C[newInd(state,i,sMP)] -= y1 * y2 * B[ind]
				C[newInd(state,i,s00)] -= y1 * x * B[ind]
			elseif sp==sMP
				C[newInd(state,i,sXX)] -= y2 * ξ * B[ind]
				C[newInd(state,i,sPM)] -= y2 * y1 * B[ind]
				C[newInd(state,i,s00)] -= y2 * x * B[ind]
			elseif sp==s00
				C[newInd(state,i,sXX)] -= x * ξ * B[ind]
				C[newInd(state,i,sPM)] -= x * y1 * B[ind]
				C[newInd(state,i,sMP)] -= x * y2 * B[ind]
			elseif sp==s0P
				C[newInd(state,i,sP0)] -= y1 * y2 * B[ind]
				C[newInd(state,i,sMM)] -= y1 * z * B[ind]
			elseif sp==sP0
				C[newInd(state,i,s0P)] -= y2 * y1 * B[ind]
				C[newInd(state,i,sMM)] -= y2 * z * B[ind]
			elseif sp==sMM
				C[newInd(state,i,s0P)] -= z * y1 * B[ind]
				C[newInd(state,i,sP0)] -= z * y2 * B[ind]
			elseif sp==s0M
				C[newInd(state,i,sM0)] -= y2 * y1 * B[ind]
				C[newInd(state,i,sPP)] -= y2 * z * B[ind]
			elseif sp==sM0
				C[newInd(state,i,s0M)] -= y1 * y2 * B[ind]
				C[newInd(state,i,sPP)] -= y1 * z * B[ind]
			elseif sp==sPP
				C[newInd(state,i,s0M)] -= z * y2 * B[ind]
				C[newInd(state,i,sM0)] -= z * y1 * B[ind]
			end
		end
	end
end

println("preparing...")
@time for i = 1 : 4^L
	setFlag!(flag_,i)
	computeDiag!(diag_,flag_,i)
end


println("computing eigenvalues...")
H=LinearMap((C,B)->Hfunc!(C,B,diag_,flag_),4^L,ismutating=true,issymmetric=true,isposdef=false)
@time e,v = eigs(H,nev=8,which=:SR)
println(sort(e))


#=

Translation (lattice shift) stuff

=#

function Tind(ind::Int64,L::Int64,right::Bool)
	if ind==2 && iseven(L)
		return 1
	end
	if ind==1 && iseven(L)
		return 2
	end
	state = (ind-1) & (4^L-1)
	if right
		return 1+(state>>2)+((state&3)<<(2*(L-1)))
	else
		return 1+(state<<2)&(4^L-1)+(state>>(2*(L-1)))
	end
end

function Tfunc!(C::Vector{Float64},B::Vector{Float64},L::Int64=L,right::Bool=true)
	for ind = 1 : 4^L
		C[ind] = B[Tind(ind,L,right)]
	end
end

T=LinearMap((C,B)->Tfunc!(C,B),4^L,ismutating=true,issymmetric=false,isposdef=false)


# Print as mathematica array to reuse Mathematica code for making plots.
function mathematicaVector(V::Vector{Float64})
	s="{"
	for i = 1:(size(V,1)-1)
		s*=string(V[i])
		s*=", "
	end
	s*=string(V[size(V,1)])
	s*="}"
	return s
end

# Print as mathematica matrix to reuse Mathematica code for making plots.
function mathematicaMatrix(M::Vector{Vector{Float64}})
	s="{\n"
	for i = 1:(size(M,1)-1)
		s*=mathematicaVector(M[i])
		s*=",\n"
	end
	s*=mathematicaVector(M[size(M,1)])
	s*="\n}\n"
	return s
end

#=
Simultaneously diagonalize Hamiltonian and translation.
Output sorted pairs of energy and momentum.
=#
println()
smallH = Matrix(diagm(e))
smallT = Matrix(adjoint(v)*T*v)
smalle,smallv = eigen(smallH+smallT)
Hs = real(diag(adjoint(smallv)*smallH*smallv))
Ts = complex(diag(adjoint(smallv)*smallT*smallv))
Ps = real(map(log,Ts)/(2*π*im))*L
HPs = hcat(Hs,Ps)
HPs = sort([HPs[i,:] for i in 1:size(HPs, 1)])
s=""
s*=string(HPs[1][1])
println(mathematicaMatrix(HPs))

#=
A fusionFlag retains the information of whether main flag is 0 to save memory.
=#
function setFusionFlag!(flag::Vector{Bool},ind::Int64,L::Int64=L+2)
	state=ind-1
	below=(state>>(2*(L+1)))
	start=(state>>(2*L))&3
	state=state&(4^L-1)
	if (start==3) || (below>=L)
		flag[ind] = true
		return
	end
	if state==0 && isodd(L)
		flag[ind] = true
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
			evenxs = ! evenxs
			if ((pos+1)==below || (below!=0 && pos==L-1))
				tot = 3-tot
			end
		else
			if(!evenxs)
				flag[ind] = true
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
		flag[ind] = true
		return
	end
	if tot!=0
		flag[ind] = true
	end
end

println("preparing fusion flags...")

fusionFlag_ = zeros(Bool,4^L)
@time for i = 1 : 4^L
	setFusionFlag!(fusionFlag_,i,L)
end

extendedFusionFlag_ = zeros(Bool,4^(L+3)*(L+2))
@time for i = 1 : 4^(L+3)*(L+2)
	setFusionFlag!(extendedFusionFlag_,i,L+2)
end

#=
Distinguished from Yuji's mainFlag by variable type.
mainFlag(flag, ...) and mainFlag(fusionFlag, ...) return the same thing.
=#
function mainFlag(flag::Vector{Bool},ind::Int64,L::Int64=L)::Bool
	return flag[ind]
end

#=

An edgeState is a 6^L encoding of a state using the original anyon labels.
An edgeStateMapping maps from index (4^L vertex label encoding) to edgeState.
This is convenient for state visualization and debugging.

It is also used by zip! to infer the correct F-symbol.
But for this purpose, it suffices to know the anyon label before the draped ρ.
TODO Only retaining this minimal info will save memory.
How many nearby anyon labels to keep (1,2,3) depends on time vs memory.

=#

const edgeMapping=Dict('1'=>1, 'a'=>2, 'b'=>3, 'ρ'=>4, 'σ'=>5, 'τ'=>6)
const edgeRevmapping=Dict(1=>'1', 2=>'a', 3=>'b', 4=>'ρ', 5=>'σ', 6=>'τ')

function setEdgeStateMapping!(edgeStateMapping::Vector{Int64},ind::Int64,flag::Vector{Bool},L::Int64=L)
	below = ((ind-1)>>(2*(L+1)))
	start = ((ind-1)>>(2*L)) & 3
	state = (ind-1)&(4^L-1)
	evenxs=iseven(trailingXs(state,L))
	if(state==1 && iseven(L))
		state=0
		evenxs=false
	end
	tot=start
	if evenxs
		edgeStateMapping[ind] |= 4+start
	else
		edgeStateMapping[ind] |= 1+start
	end
	for pos = 0 : L-2
		a=(state >> (2*pos)) & 3
		if a==0
			evenxs = ! evenxs
			if ((pos+1)==below || (below!=0 && pos==L-1))
				tot = 3-tot
			end
		elseif a==2
			tot+=1
		elseif a==3
			tot+=2
		end
		tot%=3
		if evenxs
			edgeStateMapping[ind] |= ((4+tot) << (3*(pos+1)))
		else
			edgeStateMapping[ind] |= ((1+tot) << (3*(pos+1)))
		end
	end
	# revEdgeStateMapping[edgeStateMapping[ind]] = ind
end

function stringFromEdgeState(edgeState::Int64,L::Int64=L)
	s=""
	t=edgeRevmapping[(edgeState&7)]
	for i in 1 : L
		s*=edgeRevmapping[(edgeState&7)]
		edgeState>>=3
	end
	return s*t
end

# revEdgeStateMapping_ = zeros(Int64,8^L)
println("preparing edge state mapping...")
edgeStateMapping_ = zeros(Int64,4^(L+1)*L)
@time for i = 1 : 4^(L+1)*L
	setEdgeStateMapping!(edgeStateMapping_,i,fusionFlag_,L)
end

println("preparing extended edge state mapping...")
extendedEdgeStateMapping_ = zeros(Int64,4^(L+3)*(L+2))
# extendedRevEdgeStateMapping_ = zeros(Int64,8^(L+2))
@time for i = 1 : 4^(L+3)*(L+2)
	setEdgeStateMapping!(extendedEdgeStateMapping_,i,extendedFusionFlag_,L+2)
end


#=

F-symbol stuff.

=#

function isInvertible(i::Int64)
	return i<4
end

function dual(i::Int64)
	if i==2
		return 3
	elseif i==3
		return 2
	else
		return i
	end
end

function fusion(i::Int64,j::Int64)
	ans = []
	if i<4
		if j<4
			append!(ans,1+((i+j-2)%3))
		else
			append!(ans,4+((i+j-2)%3))
		end
	else
		if j<4
			return fusion(1+((4-j)%3),i)
		else
			append!(ans,1+((3+i-j)%3))
			append!(ans,[4,5,6])
		end
	end
	return ans
end

function hasFusion(i::Int64,j::Int64,k::Int64)
	fused = fusion(i,j)
	return k in fused
end

function add(i::Int64,j::Int64)
	return 4+((i+j-1)%3)
end

function FSymbol(i::Int64,j::Int64,k::Int64,l::Int64,m::Int64,n::Int64)
	if !( hasFusion(i,j,m) && hasFusion(k,dual(l),dual(m)) && hasFusion(dual(l),i,dual(n)) && hasFusion(j,k,n) )
		return 0
	end
	if isInvertible(i) || isInvertible(j) || isInvertible(k) || isInvertible(l)
		return 1
	end
	if isInvertible(m) && isInvertible(n)
		return 1/ζ
	end
	if isInvertible(m) || isInvertible(n)
		return ξ
	end
	if i!=4
		return FSymbol(4, j, add(k,i-4), l, m, n)
	end
	if j!=4
		return FSymbol(4, 4, k, add(l,j-4), m, n)
	end
	if k!=4
		return FSymbol(4, 4, 4, add(l,4-k), m, add(n,4-k))
	end
	if m!=4
		return FSymbol(4, 4, 4, l, 4, add(n,m-4))
	end
	if i==j==k==m==4 && !(isInvertible(l)) && !(isInvertible(n))
		if l+n==8
			return x
		elseif l+n==9
			return y1
		elseif l+n==10
			return y2
		elseif l+n==11
			return z
		elseif l+n==12
			return y1
		end
	end
	error("FSymbol not found")
end


#=

Zipper stuff.

=#

# Inverse of index: x->x, 0->0, ±->∓
function inv(ind::Int64)
	if ind==2
		return 3
	elseif ind==3
		return 2
	else
		return ind
	end
end

function attachInd(ind::Int64,sp::Tuple{Int64,Int64},start::Int64,L::Int64=L+2)
	if ind==2 && iseven(L)
		state = 0
	else
		state = (ind-1)
	end
	state = state << 2

	(a,b)=sp

	i = L
	state &= ~(3<<(2*(i-1)))
	state |= (a<<(2*(i-1)))

	j = 1
	state &= ~(3<<(2*(j-1)))
	state |= (b<<(2*(j-1)))

	if (state==0) && iseven(L) && (ind==1)
		state = 1
	end

	return 1+(state+(start<<(2*L))+(1<<(2*(L+1))))
end

function attach!(C::Vector{Float64},B::Vector{Float64})
	for ind = 1 : 4^(L+3)*(L+2)
		C[ind] = 0
	end
	for ind = 1 : 4^L
		if B[ind] == 0 || mainFlag(flag_,ind,L) != 0
			continue
		end
		state = stateFromInd(ind)
		if (isodd(trailingXs(state)) || (iseven(L) && ind==2)) # start label is 1
			ni = attachInd(ind,sXX,0)
			if mainFlag(extendedFusionFlag_,ni,L+2) != 0
				error("disallowed state")
			end
			C[ni] += B[ind]
		else
			ni = attachInd(ind,sXX,0)
			if mainFlag(extendedFusionFlag_,ni,L+2) != 0
				error("disallowed state")
			end
			C[ni] += 1/ζ * B[ind]
			# ni = newInd(state,L+2,s00,0,0,L+2)
			ni = attachInd(ind,s00,0)
			if mainFlag(extendedFusionFlag_,ni,L+2) != 0
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
			# ni = newInd(state,L+2,sPM,0,1,L+2)
			ni = attachInd(ind,sPM,1)
			if mainFlag(extendedFusionFlag_,ni,L+2) != 0
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
			# ni = newInd(state,L+2,sMP,0,2,L+2)
			ni = attachInd(ind,sMP,2)
			if mainFlag(extendedFusionFlag_,ni,L+2) != 0
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
		end
	end
end

function ZipInd(ind::Int64,sp::Tuple{Int64,Int64},L::Int64=L+2)
	below = ((ind-1)>>(2*(L+1)))
	if below == 0
		error("no ρ from below")
	end
	start = ((ind-1)>>(2*L)) & 3
	i = below
	below += 1
	state = (ind-1) & (4^L-1)
	if state==1 && iseven(L)
		state = 0
	end

	(a,b)=sp

	state &= ~(3<<(2*(i-1)))
	state |= (a<<(2*(i-1)))

	j=i+1
	state &= ~(3<<(2*(j-1)))
	state |= (b<<(2*(j-1)))

	if (state==0) && (isodd(trailingXs((ind-1)&(4^L-1),L)) || (iseven(L) && ((ind-1)&(4^L-1))==1))
		state = 1
	end

	return 1+(state+(start<<(2*L))+(below<<(2*(L+1))))
end

function nextEdge(e::Int64,s::Int64,below::Bool=false)
	if e<4
		if s==0
			if below
				return 4+((4-e)%3)
			else
				return e+3
			end
		else
			return false
		end
	else
		if s==0
			if below
				return 1+((7-e)%3)
			else
				return e-3
			end
		else
			return 4+((e+s-2)%3)
		end
	end
end

function zip!(C::Vector{Float64},B::Vector{Float64},i::Int64)
	for ind = 1 : 4^(L+3)*(L+2)
		C[ind] = 0
	end
	for ind = 4^(L+3)*i+1 : 4^(L+3)*(i+1)
		if B[ind] == 0 || mainFlag(extendedFusionFlag_,ind,L+2) != 0
			continue
		end

		edgeState = extendedEdgeStateMapping_[ind]
		j = i
		e1 = (edgeState>>(3*(j-1)))&7
		j += 1
		e2 = (edgeState>>(3*(j-1)))&7
		j += 1
		e3 = (edgeState>>(3*(j-1)))&7

		state = stateFromInd(ind,L+2)
		s1,s2 = localStatePair(state,i,L+2)

		if (e2 != nextEdge(e1,s1,true) || e3 != nextEdge(e2,s2,false))
			error("inconsistent edges")
		end
		for s3 = 0 : 3
			e4 = nextEdge(e1,s3,false)
			if e4 == false
				continue
			end
			for s4 = 0 : 3
				if (e3 != nextEdge(e4,s4,true))
					continue
				end
				ni = ZipInd(ind,(s3,s4))
				if mainFlag(extendedFusionFlag_,ni,L+2)!=0 || FSymbol(4,e1,4,e3,e2,e4)==0
					continue
				end
				C[ni] += FSymbol(4,e1,4,e3,e2,e4) * B[ind]
			end
		end
	end
end

function Detach!(C::Vector{Float64},B::Vector{Float64})
	for ind = 1 : 4^L
		C[ind] = 0
	end
	for ind = 4^(L+3)*(L+1)+1 : 4^(L+3)*(L+2)
		if B[ind] == 0 || mainFlag(extendedFusionFlag_,ind,L+2) != 0
			continue
		end

		state = stateFromInd(ind,L+2)
		ni = 1+(state&(4^L-1))
		if ni==1 && ((ind-1)&(4^(L+2)-1))==1 && iseven(L)
			ni=2
		end
		if mainFlag(flag_,ni,L) != 0
			continue
		end
		s1,s2 = localStatePair(state,L+1,L+2)
		if s1==inv(s2)
			if (isodd(trailingXs(state,L+2)) || iseven(L) && ni==2) # start label is 1
				if (s1,s2)==sXX
					C[ni] += ζ * B[ind]
				end
			else
				if (s1,s2)==sXX
					C[ni] += B[ind]
				else
					C[ni] += √ζ * B[ind]
				end
			end
		end
	end
end

ρ = LinearMap((C,B)->attach!(C,B),4^(L+3)*(L+2),4^L,ismutating=true,issymmetric=false,isposdef=false)
for i = 1 : L
	global ρ = LinearMap((C,B)->zip!(C,B,i),4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false) * ρ
end
ρ = LinearMap((C,B)->Detach!(C,B),4^L,4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false) * ρ

println()
smallH = Matrix(diagm(e))
smallT = Matrix(adjoint(v)*T*v)
smallρ = Matrix(adjoint(v)*ρ*v)
smalle,smallv = eigen(smallH+smallT+smallρ)
Hs = real(diag(adjoint(smallv)*smallH*smallv))
Ts = complex(diag(adjoint(smallv)*smallT*smallv))
Ps = real(map(log,Ts)/(2*π*im))*L
ρs = real(diag(adjoint(smallv)*smallρ*smallv))
HPρs = hcat(Hs,Ps,ρs)
HPρs = sort([HPρs[i,:] for i in 1:size(HPρs, 1)])
s=""
s*=string(HPρs[1][1])
print(mathematicaMatrix(HPρs))
