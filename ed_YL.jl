using LinearAlgebra,LinearMaps
using SparseArrays
using ArnoldiMethod
using Arpack
using Profile
using Traceur

const L=6
# const N=6

#=

As Julia's array is 1-based, I use

	"ind" = 1 ... 4^L

as the array index.

"state" is what you obtain by encoding edge labels using the following mapping,
and takes the values 0 ...  4^L-1.

I have decided to change to ind = state + 1.
"state" now has two meanings, either that defined above, or with additional
information about start label and below position.

=#

const mapping=Dict('x'=>0, '0'=>1, '+'=>2, '-'=>3, 'y'=>0, 'z'=>1, 'p'=>2, 'm'=>3)
const revmapping=Dict(0=>'x', 1=>'0', 2=>'+', 3=>'-')
const startMapping=Dict('0'=>0, '1'=>1, '2'=>2)

#=

The conversion between "ind" and "state" is basically done by considering modulo 4^L.

When L is even, however, some additional care is needed,
since "xxx...x" can come with two "start labels".
My convention is to use
	ind == 4^L  for "x....x" with ρ as the start label
	ind == 1  for "x....x" with 1 as the start label

=#

function stateFromString(s::String,L=L)
	if length(s)!=L+2
		error("length doesn't match")
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

function stringFromState(state,L=L)
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

function trailingXs(state,L=L)
	state = state & (4^L-1)
	i=0
	while state!=0
		state>>=2
		i+=1
	end
	return L-i
end

function flag(ind,L=L)
	f = 0
	flagshift=L+2
	state=ind-1
	below=(state>>(2*(L+1)))
	start=(state>>(2*L))&3
	state=state&(4^L-1)
	if (start==3) || (below>=L)
		f |= 3 << flagshift
		return f
	end
	if state==0 && isodd(L)
		f |= 3 << flagshift
		return f
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
				f |= 1<<pos
			end
			evenxs = ! evenxs
			if ((pos+1)==below || (below!=0 && pos==L-1))
				tot = 3-tot
			end
		else
			if(!evenxs)
				# not allowed
				f |= 3 << flagshift
				return f
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
		f |= 3 << flagshift
		return f
	end
	f |= tot << flagshift
	return f
end

function setFlag!(flag,ind,L=L)
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
	tot=tot-start # compare with start
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

println("preparing flag...")
flag_ = zeros(Int32,4^(L+1)*L)
@time for i = 1 : 4^(L+1)*L
	setFlag!(flag_,i)
end

function setLongFlag!(flag,ind,L=L+2)
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
			# if(evenxs)
			# 	flag[ind] |= 1<<pos
			# end
			evenxs = ! evenxs
			if ((pos+1)==below || (below!=0 && pos==L-1))
				tot = 3-tot
			end
		else
			if(!evenxs)
				# not allowed
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
	tot=tot-start # compare with start
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

println("preparing new long flag...")
longFlag_ = zeros(Bool,4^(L+3)*(L+2))
for i = 1 : 4^(L+3)*(L+2)
	setLongFlag!(longFlag_,i,L+2)
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

function stateFromInd(ind,L=L)
	state=(ind-1)&(4^L-1)
	if state==1 && iseven(L)
		state=0
	end
	return state
end

function mainFlag(flag,ind,L0=L)::Int32
	if L0==L
		flagshift=L+2
		return (flag[ind] >> flagshift) & 3
	else
		return flag[ind]
	end
end

function nextSite(i,L=L)
	j=i+1
	if j==L+1
		j=1
	end
	return j
end

function localStatePair(state,i,L=L)
	j = nextSite(i,L)
	a = (state >> (2*(i-1))) & 3
	b = (state >> (2*(j-1))) & 3
	return a,b
end

function isρ1ρ(flag,ind,i)
	return ((flag[ind] >> (i-1)) & 1) == 1
end

function computeDiag!(diag,flag,ind)
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

function newInd(state,i,sp)
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

function Hfunc!(C,B,diag,flag)
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
for i = 1 : 4^L
	setFlag!(flag_,i)
	computeDiag!(diag_,flag_,i)
end

println("computing eigenvalues...")
H=LinearMap((C,B)->Hfunc!(C,B,diag_,flag_),4^L,ismutating=true,issymmetric=true,isposdef=false)
e,v = eigs(H,nev=1,which=:SR)
println(sort(e))


# Translation (lattice shift)
function Tind(ind,L,right)
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

function Tfunc!(C,B,L=L,right=true)
	for ind = 1 : 4^L
		C[ind] = B[Tind(ind,L,right)]
	end
end

T=LinearMap((C,B)->Tfunc!(C,B),4^L,ismutating=true,issymmetric=false,isposdef=false)


# Print spectrum in mathematica array to reuse mathematica code for making plots
function mathematicaVector(V)
	s="{"
	for i = 1:(size(V,1)-1)
		s*=string(V[i])
		s*=", "
	end
	s*=string(V[size(V,1)])
	s*="}"
	return s
end

function mathematicaMatrix(M)
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
Simultaneously diagonalize Hamiltonian and translation
Output sorted pairs of energy and momentum
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
print(mathematicaMatrix(HPs))


#=
Edge state stuff
=#

const edgeMapping=Dict('1'=>1, 'a'=>2, 'b'=>3, 'ρ'=>4, 'σ'=>5, 'τ'=>6)
const edgeRevmapping=Dict(1=>'1', 2=>'a', 3=>'b', 4=>'ρ', 5=>'σ', 6=>'τ')

function setEdgeState!(edgeState,ind,flag,L)
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
		edgeState[ind] |= 4+start
	else
		edgeState[ind] |= 1+start
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
			edgeState[ind] |= ((4+tot) << (3*(pos+1)))
		else
			edgeState[ind] |= ((1+tot) << (3*(pos+1)))
		end
	end
	# revEdgeState[edgeState[ind]] = ind
end

function edgeState!(ind,flag,L)
	es = 0
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
		es |= 4+start
	else
		es |= 1+start
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
			es |= ((4+tot) << (3*(pos+1)))
		else
			es |= ((1+tot) << (3*(pos+1)))
		end
	end
	return es
end

function stringFromEdgeState(edgeState,L=L)
	s=""
	t=edgeRevmapping[(edgeState&7)]
	for i in 1 : L
		s*=edgeRevmapping[(edgeState&7)]
		edgeState>>=3
	end
	return s*t
end

edgeState_ = zeros(Int64,4^(L+1)*L)
# revEdgeState_ = zeros(Int64,8^L)
println("preparing edge state...")
for i = 1 : 4^(L+1)*L
	setEdgeState!(edgeState_,i,flag_,L)
end

longEdgeState_ = zeros(Int64,4^(L+3)*(L+2))
# longRevEdgeState_ = zeros(Int64,8^(L+2))

println("preparing long edge state...")
for i = 1 : 4^(L+3)*(L+2)
	setEdgeState!(longEdgeState_,i,longFlag_,L+2)
end

# state = stateFromString("xy0xx-_2",L)
# println()
# # println(bitstring(ind))
# println(stringFromState(state,L))
# println(stringFromEdgeState(edgeState_[state+1],L))
# # println(stringFromEdgeState(EdgeState!(state+1,flag_,L),L))
# println()

# for ind = 1 : 4^(L+1)*L
# 	if (ind-1)&(4^L-1) == 1
# 		if mainFlag(flag_,ind) == 0
# 			println(stringFromState(ind-1,L), " = ", stringFromEdgeState(edgeState_[ind],L))
# 		else
# 			println(stringFromState(ind-1,L), " bad")
# 		end
# 	end
# end

#=
F-symbol stuff
=#

function isInvertible(i)
	return i<4
end

function dual(i)
	if i==2
		return 3
	elseif i==3
		return 2
	else
		return i
	end
end

function fusion(i,j)
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

function hasFusion(i,j,k)
	fused = fusion(i,j)
	return k in fused
end

function add(i,j)
	return 4+((i+j-1)%3)
end

function FSymbol(i,j,k,l,m,n)
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
	# println(i,j,k,l,m,n)
	error("FSymbol not found")
end

function inv(s)
	if s==2
		return 3
	elseif s==3
		return 2
	else
		return s
	end
end

function attachInd(ind,sp,start,L=L+2)
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

# state = stateFromString("x++x_0")
# println(stringFromState(state))
# println(stringFromState( attachInd(state+1,sMP,1)-1, L+2))


# longFlag_ = zeros(Int32,4^(L+3)*(L+2))
#
# println("preparing...")
# for i = 1 : 4^(L+3)*(L+2)
# 	setFlag!(longFlag_,i,L+2)
# end

function attach!(C,B)
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
			if mainFlag(longFlag_,ni,L+2) != 0
				println("a ", ind)
				error("disallowed state")
			end
			C[ni] += B[ind]
		else
			ni = attachInd(ind,sXX,0)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("b ",ind)
				println(flag(ni,L+2))
				println(longFlag_[ni])
				error("disallowed state")
			end
			C[ni] += 1/ζ * B[ind]
			# ni = newInd(state,L+2,s00,0,0,L+2)
			ni = attachInd(ind,s00,0)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("c ",ind)
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
			# ni = newInd(state,L+2,sPM,0,1,L+2)
			ni = attachInd(ind,sPM,1)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("d ",ind)
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
			# ni = newInd(state,L+2,sMP,0,2,L+2)
			ni = attachInd(ind,sMP,2)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("e ",ind)
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
		end
	end
end

function ZipInd(ind,sp,L=L+2)
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

function nextEdge(e,s,below=false)
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

function Zip!(C,B,i)
	for ind = 1 : 4^(L+3)*(L+2)
		C[ind] = 0
	end
	for ind = 4^(L+3)*i+1 : 4^(L+3)*(i+1)
		if B[ind] == 0 || mainFlag(longFlag_,ind,L+2) != 0
			continue
		end

		es = longEdgeState_[ind]
		j = i
		e1 = (es>>(3*(j-1)))&7
		j += 1
		e2 = (es>>(3*(j-1)))&7
		j += 1
		e3 = (es>>(3*(j-1)))&7

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
				if mainFlag(longFlag_,ni,L+2)!=0 || FSymbol(4,e1,4,e3,e2,e4)==0
					continue
				end
				C[ni] += FSymbol(4,e1,4,e3,e2,e4) * B[ind]
			end
		end
	end
end

function Detach!(C,B)
	for ind = 1 : 4^L
		C[ind] = 0
	end
	for ind = 4^(L+3)*(L+1)+1 : 4^(L+3)*(L+2)
		if B[ind] == 0 || mainFlag(longFlag_,ind,L+2) != 0
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

attach = LinearMap((C,B)->attach!(C,B),4^(L+3)*(L+2),4^L,ismutating=true,issymmetric=false,isposdef=false)
ρ = attach

for i = 1 : L
	global zip = LinearMap((C,B)->Zip!(C,B,i),4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)
	global ρ = zip * ρ
end

detach = LinearMap((C,B)->Detach!(C,B),4^L,4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)
ρ = detach * ρ

# mat = zeros(Bool,4^L)
# for i = 1 : 4^L
# 	if mainFlag(flag_,i) == 0
# 		g = zeros(Bool,4^L)
# 		g[i] = 1
# 		global mat = hcat(mat,g)
# 	end
# end
# mat = mat[:,2:end]
# ρx = adjoint(mat) * ρ * mat
# println(size(ρx))
# ex,vx = eigen(Matrix(ρx))
# println(real(ex))

# for ind = 1 : 4^L
# 	if mainFlag(flag_,ind,L)==0
# 		str = stringFromState(ind-1,L)
# 		println()
# 		println(str, " = ", stringFromEdgeState(edgeState_[ind],L))
# 		testV = zeros(4^L)
# 		testV[ind] = 1
# 		println("mainFlag: ", mainFlag(flag_,ind,L)==0)
#
# 		testU = ρ * testV
# 		for i = 1 : 4^L
# 			if testU[i] != 0
# 				if mainFlag(flag_,i,L) == 0
# 					println(stringFromState(i-1,L), " = ", stringFromEdgeState(edgeState_[i],L), " has value ", testU[i])
# 				end
# 			end
# 		end
# 	end
# end


# for ind = 1 : 4^L
# 	if mainFlag(flag_,ind,L)==0
# 		str = stringFromState(ind-1,L)
# 		println()
# 		println(str, " = ", stringFromEdgeState(edgeState_[ind],L))
# 		testV = zeros(4^L)
# 		testV[ind] = 1
# 		println("mainFlag: ", mainFlag(flag_,ind,L)==0)
#
# 		testU = attach * testV
# 		for i = 4^(L+3)+1 : 4^(L+3)*2
# 			if testU[i] != 0
# 				if mainFlag(longFlag_,i,L+2) == 0
# 					println(stringFromState(i-1,L+2), " = ", stringFromEdgeState(longEdgeState_[i],L+2), " has value ", testU[i])
# 				end
# 			end
# 		end
# 	end
# end

# for ind = 4^(L+3)*(L+1)+1 : 4^(L+3)*(L+2)
# 	if mainFlag(longFlag_,ind,L+2)==0
# 		str = stringFromState(ind-1,L+2)
# 		println()
# 		println(str, " = ", stringFromEdgeState(longEdgeState_[ind],L+2))
# 		testV = zeros(4^(L+3)*(L+2))
# 		testV[ind] = 1
# 		println("mainFlag: ", mainFlag(longFlag_,ind,L+2)==0)
#
# 		testU = detach * testV
# 		for i = 1 : 4^L
# 			if testU[i] != 0
# 				if mainFlag(flag_,i,L) == 0
# 					println(stringFromState(i-1,L), " = ", stringFromEdgeState(edgeState_[i],L), " has value ", testU[i])
# 				end
# 			end
# 		end
# 	end
# end



# zip = LinearMap((C,B)->Zip!(C,B,2),4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)
# for ind = 4^(L+3)+1 : 4^(L+3)*(L+2)
# 	# if mainFlag(longFlag_,ind,L+2)==0 && ind==1045
# 	if mainFlag(longFlag_,ind,L+2)==0 && ind==2069
# 		str = stringFromState(ind-1,L+2)
# 		println()
# 		println(str, " = ", stringFromEdgeState(longEdgeState_[ind],L+2))
# 		testV = zeros(4^(L+3)*(L+2))
# 		testV[ind] = 1
# 		println("mainFlag: ", mainFlag(longFlag_,ind,L+2)==0)
#
# 		testU = zip * testV
# 		for i = 4^(L+3)*2+1 : 4^(L+3)*(L+2)
# 			if testU[i] != 0
# 				if mainFlag(longFlag_,i,L+2) == 0
# 					println(stringFromState(i-1,L+2), " = ", stringFromEdgeState(longEdgeState_[i],L+2), " has value ", testU[i])
# 				end
# 			end
# 		end
# 	end
# end

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
