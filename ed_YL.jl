using LinearAlgebra,LinearMaps
using Arpack

const L=4
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

# function trailingXsYT(ind)
# 	i=0
# 	while ind!=0
# 		ind>>=2
# 		i+=1
# 	end
# 	L-i
# end


# flagYT_ = zeros(Int32,4^L)
#
# function setFlagYT!(flag,ind)
# 	flagshift=L+2
# 	state=ind
# 	if(ind==4^L)
# 		state=0
# 		if(isodd(L))
# 			# not allowed
# 			flag[ind] |= 3 << flagshift
# 			return
# 		end
# 	end
# 	evenxs=iseven(trailingXs(state))
# 	tot=0
# 	if(ind==1 && iseven(L))
# 		state=0
# 		evenxs=false
# 	end
# 	for pos = 0 : L-1
# 		a=(state >> (2*pos)) & 3
# 		if a==0
# 			if(evenxs)
# 				flag[ind] |= 1<<pos
# 			end
# 			evenxs = ! evenxs
# 		else
# 			if(!evenxs)
# 				# not allowed
# 				flag[ind] |= 3 << flagshift
# 				return
# 			end
# 			if a==2
# 				tot+=1
# 			elseif a==3
# 				tot+=2
# 			end
# 		end
# 	end
# 	tot%=3
# 	if ind==0 && isodd(L)
# 		flag[ind] |= 3<<flagshift
# 		return
# 	end
# 	flag[ind] |= tot << flagshift
# end

flag_ = zeros(Int32,4^(L+1)*L)

function flag!(ind,L=L)
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
	# return (f >> flagshift) & 3
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

println("preparing...")
for i = 1 : 4^(L+1)*L
	setFlag!(flag_,i)
	# setFlagYT!(flagYT_,i)
end

longFlag_ = zeros(Int32,4^(L+3)*(L+2))
for i = 1 : 4^(L+3)*(L+2)
	setFlag!(longFlag_,i,L+2)
end

println("begin")
for i = 1 : 4^L
	if i==4^L
		j = 1
	else
		j = i+1
	end
	if (flag_[j] != flag!(j))
		println(i)
	end
end
println("end")


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
	# state=(ind%(4^L))
	state=(ind-1)&(4^L-1)
	# if ind==4^L
	# 	state=0
	# end
	if state==1 && iseven(L)
		state=0
	end
	return state
end

# function stateFromIndYT(ind)
# 	state=ind
# 	if ind==4^L
# 		state=0
# 	end
# 	if ind==1 && iseven(L)
# 		state=0
# 	end
# 	return state
# end
#
# println("stateFromInd")
# for ind = 1 : 4^L
# 	if ind==4^L
# 		j=1
# 	else
# 		j=ind+1
# 	end
# 	if (stateFromInd(j) != stateFromIndYT(ind))
# 		println(ind)
# 	end
# end
# println("stateFromInd")

# Good
function mainFlag(flag,ind,L=L)::Int32
	flagshift=L+2
	return (flag[ind] >> flagshift) & 3
end

# Good
function nextSite(i,L=L)
	j=i+1
	if j==L+1
		j=1
	end
	return j
end

# Good
function localStatePair(state,i,L=L)
	j = nextSite(i,L)
	a = (state >> (2*(i-1))) & 3
	b = (state >> (2*(j-1))) & 3
	return a,b
end

# Good
function isρ1ρ(flag,ind,i)
	return ((flag[ind] >> (i-1)) & 1) == 1
end

# Good
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

# function newInd(state,i,sp,L=L,start=false,below=false)
# 	if start==false
# 		start = (state>>(2*L)) & 3
# 	end
# 	if below==false
# 		below = (state>>(2*(L+1)))
# 	end
# 	state = state & (4^L-1)
# 	evenxs = iseven(trailingXs(state,L))
#
# 	(a,b)=sp
#
# 	state &= ~(3<<(2*(i-1)))
# 	state |= (a<<(2*(i-1)))
#
# 	j=nextSite(i,L)
#
# 	state &= ~(3<<(2*(j-1)))
# 	state |= (b<<(2*(j-1)))
#
# 	if (state==0) && iseven(L) && !evenxs
# 		state = 1
# 	end
#
# 	# if(state==0)
# 	# 	if(isodd(i))
# 	# 		# state = 4^L
# 	# 	else
# 	# 		state = 1
# 	# 	end
# 	# end
#
# 	return 1+(state+(start<<(2*L))+(below<<(2*(L+1))))
# end
#
# for ind = 1 : 4^L
# 	if mainFlag(flag_, ind) == 0
# 		state = ind-1
# 		for i = 1 : 1
# 			if newInd(state,i,sXX) != newIndYT(state,i,sXX)
# 				if mainFlag(flag_, newIndYT(state,i,sXX)) == 0
# 					println("ind,i=", ind, ",", i)
# 					println( newInd(state,i,sXX), " ", newIndYT(state,i,sXX) )
# 				end
# 			end
# 		end
# 	end
# end

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
@time e,v = eigs(H,nev=8,which=:SR)
println(sort(e))


# Translation (lattice shift)
function TInd!(ind,L,right)
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

# # Translation (lattice shift)
# function TIndYL!(ind,L,right)
# 	state = ind-1
# 	if (ind==1 && iseven(L))
# 		return 4^L
# 	end
# 	if (ind==4^L && iseven(L))
# 		return 1
# 	end
# 	if right
# 		return (ind>>2)+((ind&3)<<(2*(L-1)))
# 	else
# 		return (ind<<2)&(4^L-1)+(ind>>(2*(L-1)))
# 	end
# end

function Tfunc!(C,B,L=L,right=true)
	for ind = 1 : 4^L
		C[ind] = B[TInd!(ind,L,right)]
	end
end

T=LinearMap((C,B)->Tfunc!(C,B),4^L,ismutating=true,issymmetric=false,isposdef=false)


# Print spectrum in Mathematica array to reuse Mathematica code for making plots
function MathematicaVector(V)
	s="{"
	for i = 1:(size(V,1)-1)
		s*=string(V[i])
		s*=", "
	end
	s*=string(V[size(V,1)])
	s*="}"
	return s
end

function MathematicaMatrix(M)
	s="{\n"
	for i = 1:(size(M,1)-1)
		s*=MathematicaVector(M[i])
		s*=",\n"
	end
	s*=MathematicaVector(M[size(M,1)])
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
Ts = diag(adjoint(smallv)*smallT*smallv)
Ps = real(map(log,Ts)/(2*π*im))*L
HPs = hcat(Hs,Ps)
HPs = sort([HPs[i,:] for i in 1:size(HPs, 1)])
s=""
s*=string(HPs[1][1])
print(MathematicaMatrix(HPs))


#=
Edge state stuff
=#

const edgeMapping=Dict('1'=>1, 'a'=>2, 'b'=>3, 'ρ'=>4, 'σ'=>5, 'τ'=>6)
const edgeRevmapping=Dict(1=>'1', 2=>'a', 3=>'b', 4=>'ρ', 5=>'σ', 6=>'τ')

function setEdgeState!(edgeState,revEdgeState,ind,flag,L)
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
	revEdgeState[edgeState[ind]] = ind
end


function EdgeState!(ind,flag,L)
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
	for i in 1 : L
		s*=edgeRevmapping[(edgeState&7)]
		edgeState>>=3
	end
	return s
end


edgeState_ = zeros(Int64,4^(L+1)*L)
revEdgeState_ = zeros(Int64,8^L)

println("preparing...")
for i = 1 : 4^(L+1)*L
	setEdgeState!(edgeState_,revEdgeState_,i,flag_,L)
end

longEdgeState_ = zeros(Int64,4^(L+3)*(L+2))
longRevEdgeState_ = zeros(Int64,8^(L+2))

println("preparing...")
for i = 1 : 4^(L+3)*(L+2)
	setEdgeState!(longEdgeState_,longRevEdgeState_,i,longFlag_,L+2)
end

# state = stateFromString("xy0xx-_2",L)
# println()
# # println(bitstring(ind))
# println(stringFromState(state,L))
# println(stringFromEdgeState(edgeState_[state+1],L))
# # println(stringFromEdgeState(EdgeState!(state+1,flag_,L),L))
# println()


#=
F-symbol stuff
=#

function IsInvertible!(i)
	return i<4
end

function Dual!(i)
	if i==2
		return 3
	elseif i==3
		return 2
	else
		return i
	end
end

function Fusion!(i,j)
	ans = []
	if i<4
		if j<4
			append!(ans,1+((i+j-2)%3))
		else
			append!(ans,4+((i+j-2)%3))
		end
	else
		if j<4
			return Fusion!(1+((4-j)%3),i)
		else
			append!(ans,1+((3+i-j)%3))
			append!(ans,[4,5,6])
		end
	end
	return ans
end

function HasFusion!(i,j,k)
	fused = Fusion!(i,j)
	return k in fused
end

function Add!(i,j)
	return 4+((i+j-1)%3)
end

function FSymbol!(i,j,k,l,m,n)
	if !( HasFusion!(i,j,m) && HasFusion!(k,Dual!(l),Dual!(m)) && HasFusion!(Dual!(l),i,Dual!(n)) && HasFusion!(j,k,n) )
		return 0
	end
	if IsInvertible!(i) || IsInvertible!(j) || IsInvertible!(k) || IsInvertible!(l)
		return 1
	end
	if IsInvertible!(m) && IsInvertible!(n)
		return 1/ζ
	end
	if IsInvertible!(m) || IsInvertible!(n)
		return ξ
	end
	if i!=4
		return FSymbol!(4, j, Add!(k,i-4), l, m, n)
	end
	if j!=4
		return FSymbol!(4, 4, k, Add!(l,j-4), m, n)
	end
	if k!=4
		return FSymbol!(4, 4, 4, Add!(l,4-k), m, Add!(n,4-k))
	end
	if m!=4
		return FSymbol!(4, 4, 4, l, 4, Add!(n,m-4))
	end
	if i==j==k==m==4 && !(IsInvertible!(l)) && !(IsInvertible!(n))
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


#=
Z3 stuff
=#

# const revZ3Mapping=Dict(0=>'0', 1=>'+', 2=>'-', 3=>'z')
#
# function Z3Flag!(ind,flag,L)
# 	z3 = 0
# 	if mainFlag(flag,ind,L) != 0
# 		return z3
# 	end
# 	below = ((ind-1)>>(2*(L+1)))
# 	start = ((ind-1)>>(2*L)) & 3
# 	state = (ind-1)&(4^L-1)
# 	# if(ind==4^L)
# 	# 	state=0
# 	# end
# 	evenxs = iseven(trailingXs(state,L))
# 	tot = start
# 	if(state==1 && iseven(L))
# 		state=0
# 		evenxs=false
# 	end
# 	for pos = 0 : L-1
# 		tot%=3
# 		# if evenxs
# 			z3 |= (tot << (2*pos))
# 		# else
# 			# edge to the left is invertible
# 			# z3 |= (3 << (2*pos))
# 		# end
# 		a=(state >> (2*pos)) & 3
# 		if a==0
# 			evenxs = ! evenxs
# 			if ((pos+1)==below || (below!=0 && pos==L-1))
# 				tot = 3-tot
# 			end
# 		elseif a==2
# 			tot+=1
# 		elseif a==3
# 			tot+=2
# 		end
# 	end
# 	return z3
# end
#
# function setZ3Flag!(z3Flag,ind,flag,L)
# 	if mainFlag(flag,ind,L) != 0
# 		return
# 	end
# 	below = ((ind-1)>>(2*(L+1)))
# 	start = ((ind-1)>>(2*L)) & 3
# 	state = (ind-1)&(4^L-1)
# 	# if(ind==4^L)
# 	# 	state=0
# 	# end
# 	evenxs = iseven(trailingXs(state,L))
# 	tot = start
# 	if(state==1 && iseven(L))
# 		state=0
# 		evenxs=false
# 	end
# 	for pos = 0 : L-1
# 		tot%=3
# 		if evenxs
# 			z3Flag[ind] |= (tot << (2*pos))
# 		else
# 			# edge to the left is invertible
# 			z3Flag[ind] |= (3 << (2*pos))
# 		end
# 		a=(state >> (2*pos)) & 3
# 		if a==0
# 			evenxs = ! evenxs
# 			if ((pos+1)==below || (below!=0 && pos==L-1))
# 				tot = 3-tot
# 			end
# 		elseif a==2
# 			tot+=1
# 		elseif a==3
# 			tot+=2
# 		end
# 	end
# end
#
# function stringFromZ3(z3Labels,L)
# 	s=""
# 	for i in 1 : L
# 		s*=revZ3Mapping[(z3Labels&3)]
# 		z3Labels>>=2
# 	end
# 	return s
# end
#
# flag_ = zeros(Int32,4^(L+1)*L)
# z3Flag_ = zeros(Int64,4^(L+1)*L)
#
# println("preparing...")
# for i = 1 : 4^(L+1)*L
# 	setFlag!(flag_,i)
# end
#
# println("preparing...")
# for i = 1 : 4^(L+1)*L
# 	setZ3Flag!(z3Flag_,i,flag_,L)
# end

# state = stateFromString("xy-xx0_1",L)
# println()
# println(stringFromState(state,L))
# println(stringFromZ3(Z3Flag!(state+1,flag_,L),L))
# println(stringFromEdgeState(EdgeState!(state+1,flag_,L),L))
# println(mainFlag(flag_,state+1))
# # println()
# # println(bitstring(flag!(state+1)))
# println()

# for i = 1 : 4^L
# 	if mainFlag(flag_,i,L) == 0
# 		println("state: ", stringFromIndex(i,L))
# 		println("edge:  ", stringFromEdgeState(edgeState_[i],L))
# 		# println("z3:    ", stringFromZ3(z3Flag_[i],L))
# 		# println("trail: ", isodd(trailingXs(stateFromInd(i,L),L)))
# 		println("state: ", stringFromIndex(newInd(stateFromInd(i,L) << 2,L+2,sXX,L+2),L+2))
# 		# println("edge: ", stringFromEdgeState(newInd(stateFromInd(i,L),L+2,sXX,L+2),L+2))
# 		println()
# 	end
# end

#
#
# #=
# Zipper stuff
# =#
#
# longFlag_ = zeros(Int32,4^(L+4)*L)
# # longZ3Flag_ = zeros(Int64,4^(L+2))
# longEdgeState_ = zeros(Int64,4^(L+4)*L)
# longRevEdgeState_ = zeros(Int64,8^(L+2))
#
#
# println("preparing...")
# for i = 1 : 4^(L+4)*L
# 	setFlag!(longFlag_,i,L+2)
# end
# # for i = 1 : 4^(L+2)
# # 	setZ3Flag!(longZ3Flag_,i,longFlag_,L+2)
# # end
# for i = 1 : 4^(L+4)*L
# 	setEdgeState!(longEdgeState_,longRevEdgeState_,i,longFlag_,L+2)
# end
#
#
# # for i = 1 : 4^(L+2)
# # 	if mainFlag(longFlag_,i,L+2) == 0
# # 		println("state: ", stringFromIndex(i,L+2))
# # 		println("edge:  ", stringFromEdgeState(longEdgeState_[i],L+2))
# # 		# println("z3:    ", stringFromZ3(z3Flag_[i],L))
# # 		println()
# # 	end
# # end
#
#
# const ζinvLabels = [(0, 0, 0, 0), (0, 0, 1, 0)]
#
# const ξLabels = [(0, 0, 0, 1), (0, 0, 0, 2), (0, 0, 0, 3), (0, 0, 1, 1), (0, 0, 1, 2),
# (0, 0, 1, 3), (0, 1, 1, 0), (0, 2, 3, 0), (0, 3, 2, 0)]
#
# const xLabels = [(0, 1, 1, 1), (0, 2, 3, 3), (0, 3, 2, 2), (1, 1, 2, 2), (1, 2, 1, 1),
# (1, 3, 3, 3), (2, 1, 3, 3), (2, 2, 2, 2), (2, 3, 1, 1)]
#
# const zLabels = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 2), (0, 2, 2, 1), (0, 3, 1, 3),
# (0, 3, 3, 1), (1, 1, 1, 3), (1, 1, 3, 1), (1, 2, 2, 3), (1, 2, 3, 2),
# (1, 3, 1, 2), (1, 3, 2, 1), (2, 1, 1, 2), (2, 1, 2, 1), (2, 2, 1, 3),
# (2, 2, 3, 1), (2, 3, 2, 3), (2, 3, 3, 2)]
#
# const y1Labels = [(0, 1, 1, 2), (0, 1, 2, 1), (0, 1, 3, 3), (0, 2, 1, 3), (0, 2, 2, 2),
# (0, 2, 3, 1), (0, 3, 1, 1), (0, 3, 2, 3), (0, 3, 3, 2), (1, 1, 1, 1),
# (1, 1, 2, 3), (1, 1, 3, 2), (1, 2, 1, 2), (1, 2, 2, 1), (1, 2, 3, 3),
# (1, 3, 1, 3), (1, 3, 2, 2), (1, 3, 3, 1), (2, 1, 1, 3), (2, 1, 2, 2),
# (2, 1, 3, 1), (2, 2, 1, 1), (2, 2, 2, 3), (2, 2, 3, 2), (2, 3, 1, 2),
# (2, 3, 2, 1), (2, 3, 3, 3)]
#
# const y2Labels = [(0, 1, 1, 3), (0, 1, 2, 2), (0, 1, 3, 1), (0, 2, 1, 1), (0, 2, 2, 3),
# (0, 2, 3, 2), (0, 3, 1, 2), (0, 3, 2, 1), (0, 3, 3, 3), (1, 1, 1, 2),
# (1, 1, 2, 1), (1, 1, 3, 3), (1, 2, 1, 3), (1, 2, 2, 2), (1, 2, 3, 1),
# (1, 3, 1, 1), (1, 3, 2, 3), (1, 3, 3, 2), (2, 1, 1, 1), (2, 1, 2, 3),
# (2, 1, 3, 2), (2, 2, 1, 2), (2, 2, 2, 1), (2, 2, 3, 3), (2, 3, 1, 3),
# (2, 3, 2, 2), (2, 3, 3, 1)]
#
function Inv!(s)
	if s==2
		return 3
	elseif s==3
		return 2
	else
		return s
	end
end
#
# Twisted Z3 charge
# function Q!(s)
# 	if s==0
# 		return 0
# 	else
# 		return s-1
# 	end
# end
#

function AttachInd(ind,sp,start,L=L+2)
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

	if (state==0) && iseven(L) && ind==1
		state = 1
	end

	return 1+(state+(start<<(2*L))+(1<<(2*(L+1))))
end

# state = stateFromString("x++x_0")
# println(stringFromState(state))
# println(stringFromState( AttachInd(state+1,sMP,1)-1, L+2))


longFlag_ = zeros(Int32,4^(L+3)*(L+2))

println("preparing...")
for i = 1 : 4^(L+3)*(L+2)
	setFlag!(longFlag_,i,L+2)
end

function Attach!(C,B)
	for ind = 1 : 4^(L+3)*(L+2)
		C[ind] = 0
	end
	for ind = 1 : 4^L
		if B[ind] == 0 || mainFlag(flag_,ind,L) != 0
			continue
		end
		state = stateFromInd(ind)
		if (isodd(trailingXs(state))) # start label is 1
			ni = AttachInd(ind,sXX,0)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("a ", ind)
				error("disallowed state")
			end
			C[ni] += B[ind]
		else
			ni = AttachInd(ind,sXX,0)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("b ",ind)
				println(flag!(ni,L+2))
				println(longFlag_[ni])
				error("disallowed state")
			end
			C[ni] += 1/ζ * B[ind]
			# ni = newInd(state,L+2,s00,0,0,L+2)
			ni = AttachInd(ind,s00,0)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("c ",ind)
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
			# ni = newInd(state,L+2,sPM,0,1,L+2)
			ni = AttachInd(ind,sPM,1)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("d ",ind)
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
			# ni = newInd(state,L+2,sMP,0,2,L+2)
			ni = AttachInd(ind,sMP,2)
			if mainFlag(longFlag_,ni,L+2) != 0
				println("e ",ind)
				error("disallowed state")
			end
			C[ni] += ξ * B[ind]
		end
	end
end

#
# function ZipF!(z3, s1, s2, s3, s4)
# 	if z3==3
# 		# first edge is invertible
# 		if s1==s3==x && s2==s4
# 			return 1
# 		else
# 			return 0
# 		end
# 	end
#
# 	if ((Q!(s1)+Q!(s2)-Q!(s3)-Q!(s4))%3) != 0 || ((s1!=0)+(s2!=0)+(s3!=0)+(s4!=0)) != 0
# 		# s1,s2 and s3,s4 must get to the same third edge
# 		return 0
# 	# elseif s1==s2==0
# 	# 	if s3!=Inv!(s4)
# 	# 		return 0
# 	# 	end
# 	# elseif (s1==0 || s2==0) && (s1!=s3 || s2!=s4)
# 	# 	return 0
# 	# elseif s3==s4==0
# 	# 	if s1!=Inv!(s2)
# 	# 		return 0
# 	# 	end
# 	# elseif (s3==0 || s4==0) && (s1!=s3 || s2!=s4)
# 	# 	return 0
# 	end
#
# 	label = (z3, s1, s2, s3)
# 	if label in ζinvLabels
# 		return 1/ζ
# 	elseif label in ξLabels
# 		return ξ
# 	elseif label in xLabels
# 		return x
# 	elseif label in zLabels
# 		return z
# 	elseif label in y1Labels
# 		return y1
# 	elseif label in y2Labels
# 		return y2
# 	else
# 		return 1
# 	end
# end
#

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

	if (state==0) && iseven(L) && ((ind-1)&(4^L-1))==1
		state = 1
	end

	return 1+(state+(start<<(2*L))+(below<<(2*(L+1))))
end

# state = stateFromString("yxxx-_1",L+2)
# println(mainFlag(longFlag_,state+1,L+2)==0)
# println(stringFromState(state,L+2), " = ", stringFromEdgeState(EdgeState!(state+1,longFlag_,L+2),L+2))
# println(stringFromState(ZipInd(state+1,sPM)-1,L+2), " = ", stringFromEdgeState(EdgeState!(ZipInd(state+1,sPM),longFlag_,L+2),L+2))

function nextEdge!(e,s,below=false)
	if e<4
		if s==0
			if below
				return 4+((4-e)%3)
			else
				return e+3
			end
		else
			# return 1+((e+s-2)%3)
			# println(e,s,below)
			# error("not allowed")
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
	for ind = 1 : 4^(L+3)*(L+2)
		if B[ind] == 0 || mainFlag(longFlag_,ind,L+2) != 0
			continue
		end
		# z3 = (longZ3Flag[ind] >> 2*(i-1)) & 3

		es = longEdgeState_[ind]
		j = i
		e1 = (es>>(3*(j-1)))&7
		j += 1
		e2 = (es>>(3*(j-1)))&7
		j += 1
		e3 = (es>>(3*(j-1)))&7

		# e1 = edgeMapping[es[i]]
		# e2 = edgeMapping[es[i+1]]
		# e3 = edgeMapping[es[i+2]]

		# if (ind == 4^L)
		# 	println(i)
		# 	println(stringFromIndex(ind,L+2))
		# 	println(localStatePair(state,3,L+2))
		# 	println()
		# end

		state = stateFromInd(ind,L+2)
		s1,s2 = localStatePair(state,i,L+2)

		# println(s1,s2)
		# println(stringFromState(ind-1,L+2), " = ", stringFromEdgeState(longEdgeState_[ind],L+2))
		if (e2 != nextEdge!(e1,s1,true) || e3 != nextEdge!(e2,s2,false))
			error("inconsistent edges")
		end
		for s3 = 0 : 3
			e4 = nextEdge!(e1,s3,false)
			if e4 == false
				continue
			end
			for s4 = 0 : 3
				if (e3 != nextEdge!(e4,s4,true))
					# error("haha")
					continue
				end
				# if e1<4
				# 	# if s3 != 0
				# 	# 	error("here")
				# 	# end
				# 	if s3 != 0
				# 		continue
				# 	end
				# 	e4=e1+3
				# else
				# 	if s3==0
				# 		e4=e1-3
				# 	else
				# 		e4=4+((e1+s3-2)%3)
				# 	end
				# end
				# if e4 != e4x
				# 	error("xxx")
				# end
				ni = ZipInd(ind,(s3,s4))
				if mainFlag(longFlag_,ni,L+2)!=0 || FSymbol!(4,e1,4,e3,e2,e4)==0
					continue
				end
				# println(4,e1,4,e3,e2,e4," = ",FSymbol!(4,e1,4,e3,e2,e4))
				C[ni] += FSymbol!(4,e1,4,e3,e2,e4) * B[ind]
				# C[ni] += ZipF!(z3,s1,s2,s3,s4) * B[ind]
			# # if (ind == 4^L)
			# 	edgeState = longEdgeState_[ind]
			# 	j = i
			# 	e1 = (edgeState>>(3*(j-1)))&7
			# 	j = nextSite(j,L+2)
			# 	e2 = (edgeState>>(3*(j-1)))&7
			# 	j = nextSite(j,L+2)
			# 	e3 = (edgeState>>(3*(j-1)))&7
			# 	if e1<4
			# 		e4=e1+3
			# 	else
			# 		if s3==0
			# 			e4=e1-3
			# 		else
			# 			e4=4+((e1+s3-2)%3)
			# 		end
			# 	end
			# 	if (ZipF!(z3,s1,s2,s3,s4) != 0) && (FSymbol!(4,e1,4,e3,e2,e4) != ZipF!(z3,s1,s2,s3,s4))
			# 		# println(i)
			# 		println(stringFromIndex(ind,L+2)," at position ",i)
			# 		# println(revmapping[s4])
			# 		println("e1...e4: ",edgeRevmapping[e1],edgeRevmapping[e2],edgeRevmapping[e3],edgeRevmapping[e4])
			# 		# println(4,e1,4,e3,e2,e4)
			# 		println("FSymbol(",4,e1,4,e3,e2,e4,") = ", FSymbol!(4,e1,4,e3,e2,e4))
			# 		# println(edgeState[nextSite(i,L+2)])
			# 		println(stringFromEdgeState(edgeState,L+2))
			# 		# println(z3,s1,s2,s3,s4)
			# 		println("ZipF(",z3,s1,s2,s3,s4,") = ",ZipF!(z3,s1,s2,s3,s4))
			# 		println()
			# 	end
				# end

				# if B[ind] != 0
				# 	println(newInd(state,i,(s3,s4),L+2))
				# end
				# if B[ind] != 0 && ZipF!(z3,s1,s2,s3,s4) != 0
				# 	println("z3,s1,s2,s3,s4: ", z3,s1,s2,s3,s4)
				# 	println("F: ", ZipF!(z3,s1,s2,s3,s4))
				# 	ni = newInd(state,i,(s3,s4),L+2)
				# 	println(stringFromIndex(ni,L+2))
				# 	println(stringFromEdgeState(longEdgeState_[ni],L+2))
				# 	println("value: ", ZipF!(z3,s1,s2,s3,s4) * B[ind])
				# 	C[ni] += ZipF!(z3,s1,s2,s3,s4) * B[ind]
				# 	println(norm(C))
				# end
			end
		end
	end
	# println(C[newInd(state,i,(0,0),L+2)])
	# println("norm: ", norm(C))
end

function Detach!(C,B)
	for ind = 1 : 4^L
		C[ind] = 0
	end
	for ind = 1 : 4^(L+3)*(L+2)
		if B[ind] == 0 || mainFlag(longFlag_,ind,L+2) != 0
			continue
		end

		# below = ((ind-1)>>(2*(L+3)))
		# if below != L+1
		# 	error("bad")
		# end

		state = stateFromInd(ind,L+2)
		ni = 1+(state&(4^L-1))
		if ni==1 && ((ind-1)&(4^(L+2)-1))==1 && iseven(L)
			ni=2
		end
		if mainFlag(flag_,ni,L) != 0
			continue
		end
		s1,s2 = localStatePair(state,L+1,L+2)
		if s1==Inv!(s2)
			if (isodd(trailingXs(state,L+2))) # start label is 1
				if (s1,s2)==sXX
					# println(stringFromIndex(state,L+2))
					# println(stringFromIndex(state&(4^L-1),L))
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

# state = stateFromString("z0xx0_1",L+2)
# println()
# println(mainFlag(longFlag_,state+1,L+2)==0)
# println(stringFromState(state,L+2), " = ", stringFromEdgeState(EdgeState!(state+1,longFlag_,L+2),L+2))
# B = zeros(4^(L+3)*(L+2))
# B[state+1] = 1
# C = zeros(4^(L+3)*(L+2))
# Zip!(C,B,1)
# println(norm(C))


detach = LinearMap((C,B)->Detach!(C,B),4^L,4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)

# state = stateFromString("00xx0_0",L+2)
# println(mainFlag(longFlag_,state+1,L+2)==0)
# println(stringFromState(state,L+2), " = ", stringFromEdgeState(EdgeState!(state+1,longFlag_,L+2),L+2))
# testV = zeros(4^(L+3)*(L+2))
# testV[state+1] = 1
# testU = detach * testV
# for i = 1 : 4^L
# 	if testU[i] != 0
# 		if mainFlag(flag_,i,L) == 0
# 			println(stringFromState(i-1), " = ", stringFromEdgeState(edgeState_[i],L), " has value ", testU[i])
# 		else
# 			println("bad flag: ", i, " = ", stringFromIndex(i))
# 		end
# 	end
# end

zip = LinearMap((C,B)->Zip!(C,B,1),4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)

# state = stateFromString("z0xx0_0",L+2)
# if mainFlag(longFlag_,state+1,L+2)==0
# 	println()
# 	println(stringFromState(state,L+2), " = ", stringFromEdgeState(EdgeState!(state+1,longFlag_,L+2),L+2))
# 	testV = zeros(4^(L+3)*(L+2))
# 	testV[state+1] = 1
# 	testU = zip * testV
# 	for i = 1 : 4^(L+3)*(L+2)
# 		if testU[i] != 0
# 			if mainFlag(longFlag_,i,L+2) == 0
# 				println(stringFromState(i-1,L+2), " = ", stringFromEdgeState(longEdgeState_[i],L+2), " has value ", testU[i])
# 			else
# 				println("bad flag: ", i, " = ", stringFromIndex(i,L+2))
# 			end
# 		end
# 	end
# end

# println(stringFromState( Ind(state+1,sMP,1)-1, L+2))

#
attach = LinearMap((C,B)->Attach!(C,B),4^(L+3)*(L+2),4^L,ismutating=true,issymmetric=false,isposdef=false)
ρ = attach

# #
for i = 1 : L
	global zip = LinearMap((C,B)->Zip!(C,B,i),4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)
	global ρ = zip * ρ
end
#
# # t = LinearMap((C,B)->Tfunc!(C,B,L+2,false),4^(L+2),ismutating=true,issymmetric=false,isposdef=false)
# # ρ = t * ρ
#

detach = LinearMap((C,B)->Detach!(C,B),4^L,4^(L+3)*(L+2),ismutating=true,issymmetric=false,isposdef=false)
ρ = detach * ρ

# println("computing eigenvalues...")
# @time e,v = eigs(ρ,nev=1,which=:SR)
# println(sort(real(e)))

# println(size(ρ))
#
# # println(Matrix(adjoint(v) * ρ * v))
#
# println(norm(Matrix(ρ)))
# good = zeros(Bool,4)
mat = zeros(Bool,4^L)
for i = 1 : 4^L
	if mainFlag(flag_,i) == 0
		g = zeros(Bool,4^L)
		g[i] = 1
		global mat = hcat(mat,g)
	end
end
mat = mat[:,2:end]

ρ = adjoint(mat) * ρ * mat
println(size(ρ))
ex,vx = eigen(Matrix(ρ))
println(real(ex))

#
#
# # edgeState = longEdgeState_
# # flag = longFlag_
#
#
# # str = "000"
# # ind = indexFromString(str,inL)
#
# for ind = 1 : 4^(L+3)*(L+2)
#
# 	inL = L+2
# 	edgeState = longEdgeState_
# 	flag = longFlag_
#
# 	if mainFlag(flag,ind,inL)==0
#
# 		str = stringFromState(ind-1,inL)
# 		println()
# 		println(str, " = ", stringFromEdgeState(edgeState[ind],inL))
# 		testV = zeros(4^(L+3)*(L+2))
# 		testV[ind] = 1
# 		println("mainFlag: ", mainFlag(flag,ind,inL)==0)
#
# 		# global zip = LinearMap((C,B)->Zip!(C,B,1,longZ3Flag_),4^(L+2),ismutating=true,issymmetric=false,isposdef=false)
# 		# testU = zip * testV
#
# 		testU = detach * testV
#
# 		outL = L
# 		flag = flag_
# 		# edgeState = edgeState_
# 		edgeState = edgeState_
# 		for i = 1 : 4^L
# 			if testU[i] != 0
# 				if mainFlag(flag,i,outL) == 0
# 					# println(edgeState[ind])
# 					# println(stringFromState(i-1,outL), " has value ", testU[i])
# 					println(stringFromState(i-1,outL), " = ", stringFromEdgeState(edgeState[i],outL), " has value ", testU[i])
# 				# # else
# 				# 	println("bad flag: ", ind, " = ", stringFromIndex(ind,outL))
# 				end
# 			end
# 		end
#
# 	end
# end
#
#

# println()
# smallH = Matrix(diagm(e))
# smallT = Matrix(adjoint(v)*T*v)
# smallρ = Matrix(adjoint(v)*ρ*v)
# smalle,smallv = eigen(smallH+smallT+smallρ)
# Hs = real(diag(adjoint(smallv)*smallH*smallv))
# Ts = diag(adjoint(smallv)*smallT*smallv)
# Ps = real(map(log,Ts)/(2*π*im))*L
# ρs = real(diag(adjoint(smallv)*smallρ*smallv))
# HPρs = hcat(Hs,Ps,ρs)
# HPρs = sort([HPρs[i,:] for i in 1:size(HPρs, 1)])
# s=""
# s*=string(HPρs[1][1])
# print(MathematicaMatrix(HPρs))




# HPs = sort([HPs[i,:] for i in 1:size(HPs, 1)])
# s=""
# s*=string(HPs[1][1])
# print(MathematicaMatrix(HPs))

#
#
# # function Zip!(C,B,flag,z3Flag)
# # 	for ind = 1 : 4^(L+2)
# # 		state = stateFromIndLong(ind)
# # 		for i = 1 : L
# # 			sp=localStatePair(state,i)
# # 			if sp==sXX  && isρ1ρ(flag,ind,i)
# # 				C[newIndLong(state,i,sPM)] -= ξ * y1 * B[ind]
# # 				C[newIndLong(state,i,sMP)] -= ξ * y2 * B[ind]
# # 				C[newIndLong(state,i,s00)] -= ξ * x * B[ind]
# # 			elseif sp==sPM
# # 				C[newIndLong(state,i,sXX)] -= y1 * ξ * B[ind]
# # 				C[newIndLong(state,i,sMP)] -= y1 * y2 * B[ind]
# # 				C[newIndLong(state,i,s00)] -= y1 * x * B[ind]
# # 			elseif sp==sMP
# # 				C[newIndLong(state,i,sXX)] -= y2 * ξ * B[ind]
# # 				C[newIndLong(state,i,sPM)] -= y2 * y1 * B[ind]
# # 				C[newIndLong(state,i,s00)] -= y2 * x * B[ind]
# # 			elseif sp==s00
# # 				C[newIndLong(state,i,sXX)] -= x * ξ * B[ind]
# # 				C[newIndLong(state,i,sPM)] -= x * y1 * B[ind]
# # 				C[newIndLong(state,i,sMP)] -= x * y2 * B[ind]
# # 			elseif sp==s0P
# # 				C[newIndLong(state,i,sP0)] -= y1 * y2 * B[ind]
# # 				C[newIndLong(state,i,sMM)] -= y1 * z * B[ind]
# # 			elseif sp==sP0
# # 				C[newIndLong(state,i,s0P)] -= y2 * y1 * B[ind]
# # 				C[newIndLong(state,i,sMM)] -= y2 * z * B[ind]
# # 			elseif sp==sMM
# # 				C[newIndLong(state,i,s0P)] -= z * y1 * B[ind]
# # 				C[newIndLong(state,i,sP0)] -= z * y2 * B[ind]
# # 			elseif sp==s0M
# # 				C[newIndLong(state,i,sM0)] -= y2 * y1 * B[ind]
# # 				C[newIndLong(state,i,sPP)] -= y2 * z * B[ind]
# # 			elseif sp==sM0
# # 				C[newIndLong(state,i,s0M)] -= y1 * y2 * B[ind]
# # 				C[newIndLong(state,i,sPP)] -= y1 * z * B[ind]
# # 			elseif sp==sPP
# # 				C[newIndLong(state,i,s0M)] -= z * y2 * B[ind]
# # 				C[newIndLong(state,i,sM0)] -= z * y1 * B[ind]
# # 			else
# # 				return
# # 			end
# # 		end
# # 	end
# # end
