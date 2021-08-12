using LinearAlgebra,LinearMaps
using Arpack

const L=3
# const N=6

#=

As Julia's array is 1-based, I use

	"ind" = 1 ... 4^L

as the array index.

"state" is what you obtain by encoding edge labels using the following mapping,
and takes the values 0 ...  4^L-1.

=#

const mapping=Dict('x'=>0, '0'=>1, '+'=>2, '-'=>3, 'y'=>0)
const revmapping=Dict(0=>'x', 1=>'0', 2=>'+', 3=>'-')

#=

The conversion between "ind" and "state" is basically done by considering modulo 4^L.

When L is even, however, some additional care is needed,
since "xxx...x" can come with two "start labels".
My convention is to use
	ind == 4^L  for "x....x" with ρ as the start label
	ind == 1  for "x....x" with 1 as the start label

=#

function indexFromString(s::String,start=0,L=L)
	if length(s)!=L
		error("length doesn't match")
	end
	a=start<<(2*(L+1))
	# cnt=0
	for i in L : -1 : 1
		if s[i]=='y'
			a+=(i<<(2*(L+2)))
			# cnt+=1
		end
		a+=(mapping[s[i]]<<(2*(i-1)))
	end
	if (a%(4^L))==0
		a+=4^L
	end
	return a
end

function stringFromIndex(ind,L=L)
	below = (ind>>(2*(L+2)))
	start = (ind>>(2*(L+1))) & 3
	# if (ind>>(2*L))==0
	# 	below = 0
	# else
	# 	below = [ind>>(2*L+2)&(2^N-1))]
	# end
	ind=ind%(4^L)
	# if ind==4^L
	# 	ind=0
	# end
	s=string(start)
	for i in 1 : L
		if i == below
			s*='y'
		else
			s*=revmapping[(ind&3)]
		end
		ind>>=2
	end
	return s
end


function trailingXs(ind,L=L)
	ind=ind%(4^L)
	i=0
	while ind!=0
		ind>>=2
		i+=1
	end
	L-i
end

# const flagshift = L+2


function bitdump(i,L=L)
	flagshift=L+2
	s=""
	for pos = 0 : L-1
		if (i & (1<<pos))!=0
			s*="1"
		else
			s*="0"
		end
	end
	println(s)
	a = (i >> flagshift) & 3
	if(a==3)
		println("not allowed")
	else
		println("Z3 charge: $a")
	end
end


#=

flag[ind] is a bitmask; from the 0th bit to (L-1)-th bit , 1 indicates that
the edge labels (i+1) (i+2) are x,x and corresponds to ρ,1,ρ
the (L+2)th and (L+3)th bit combine to form 0,1,2,3 ,
where 0,1,2 are twisted Z3 charge and 3 means it is a forbidden state

=#

flag_ = zeros(Int32,4^L)

function setFlag!(flag,ind,L=L)
	flagshift = L+2
	below=(ind>>(2*(L+2)))
	if below != 0
		return
	end
	state=(ind%(4^L))
	if ind==0 && isodd(L)
		flag[ind] |= 3 << flagshift
		return
	end
	# if(ind==4^L)
	# 	state=0
	# 	if(isodd(L))
	# 		# not allowed
	# 		flag[ind] |= 3 << flagshift
	# 		return
	# 	end
	# end
	evenxs=iseven(trailingXs(state,L))
	tot=0
	if(ind==1 && iseven(L))
		state=0
		evenxs=false
	end
	for pos = 0 : L-1
		a=(state >> (2*pos)) & 3
		if a==0
			if(evenxs)
				flag[ind] |= 1<<pos
			end
			evenxs = ! evenxs
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
	tot%=3
	if ind==0 && isodd(L)
		flag[ind] |= 3<<flagshift
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

function stateFromInd(ind,L=L)
	state=(ind%(4^L))
	if ind==4^L
		state=0
	end
	if ind==1 && iseven(L)
		state=0
	end
	return state
end

function mainFlag(flag,ind,L=L)::Int32
	flagshift=L+2
	return (flag[ind] >>flagshift) & 3
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



function newInd(state,i,sp,L=L)
	(a,b)=sp

	state &= ~(3<<(2*(i-1)))
	state |= (a<<(2*(i-1)))

	j=nextSite(i,L)

	state &= ~(3<<(2*(j-1)))
	state |= (b<<(2*(j-1)))

	if(state!=0)
		return state
	end
	if(isodd(i))
		return 4^L
	else
		return 1
	end
end

function pettyPrint(v)
	for x in v
		if(abs(x)<.0001)
			x=0
		end
		print("$x,")
	end
	println("")
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


# Translation (lattice shift)
function TInd!(ind,L,right)
	if (ind==1 && iseven(L))
		return 4^L
	end
	if (ind==4^L && iseven(L))
		return 1
	end
	if right
		return (ind>>2)+((ind&3)<<(2*(L-1)))
	else
		return (ind<<2)&(4^L-1)+(ind>>(2*(L-1)))
	end
end

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


# println("computing eigenvalues...")
# H=LinearMap((C,B)->Hfunc!(C,B,diag_,flag_),4^L,ismutating=true,issymmetric=true,isposdef=false)
# @time e,v = eigs(H,nev=8,which=:SR)
# println(sort(e))
#
#=
Simultaneously diagonalize Hamiltonian and translation
Output sorted pairs of energy and momentum
=#
# println()
# smallH = Matrix(diagm(e))
# smallT = Matrix(adjoint(v)*T*v)
# smalle,smallv = eigen(smallH+smallT)
# Hs = real(diag(adjoint(smallv)*smallH*smallv))
# Ts = diag(adjoint(smallv)*smallT*smallv)
# Ps = real(map(log,Ts)/(2*π*im))*L
# HPs = hcat(Hs,Ps)
# HPs = sort([HPs[i,:] for i in 1:size(HPs, 1)])
# s=""
# s*=string(HPs[1][1])
# print(MathematicaMatrix(HPs))


#=
Edge state stuff
=#
const edgeRevmapping=Dict(1=>'1', 2=>'a', 3=>'b', 4=>'ρ', 5=>'σ', 6=>'τ')

function setEdgeState!(edgeState,revEdgeState,ind,flag,L)
	below = (ind>>(2*(L+2)))
	start = (ind>>(2*(L+1))) & 3
	# println("start: ",start)
	# ind = ind%(4^L)
	state = ind%(4^L)
	# if ind==0
	# 	ind=4^L
	# end
	# if mainFlag(flag,ind,L) != 0
	# 	return
	# end
	# state=ind
	# if(ind==4^L)
	# 	state=0
	# end
	evenxs=iseven(trailingXs(state,L))
	if(ind==1 && iseven(L))
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
			if (pos+1)==below && evenxs
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
	below = (ind>>(2*(L+2)))
	start = (ind>>(2*(L+1))) & 3
	state = ind%(4^L)
	evenxs=iseven(trailingXs(state,L))
	if(ind==1 && iseven(L))
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
			if (pos+1)==below && evenxs
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


function stringFromEdgeState(edgeState,L)
	s=""
	for i in 1 : L
		s*=edgeRevmapping[(edgeState&7)]
		edgeState>>=3
	end
	return s
end


edgeState_ = zeros(Int64,4^(L+2)*L)
revEdgeState_ = zeros(Int64,8^L)

println("preparing...")
for i = 1 : 4^(L+2)*L
	setEdgeState!(edgeState_,revEdgeState_,i,flag_,L)
end

ind = indexFromString("y0x",2,L)
println()
# println(bitstring(ind))
println(stringFromIndex(ind,L))
println(stringFromEdgeState(edgeState_[ind],L))
println(stringFromEdgeState(EdgeState!(ind,flag_,L),L))
println()

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
	if i==j==k==l==4 && !(IsInvertible!(m)) && !(IsInvertible!(n))
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
Z3 stuff
=#

const revZ3Mapping=Dict(0=>'0', 1=>'+', 2=>'-', 3=>'1')

z3Flag_ = zeros(Int64,4^L)

function setZ3Flag!(z3Flag,ind,flag,L)
	if mainFlag(flag,ind,L) != 0
		return
	end
	state=ind
	if(ind==4^L)
		state=0
	end
	evenxs=iseven(trailingXs(state,L))
	tot=0
	if(ind==1 && iseven(L))
		state=0
		evenxs=false
	end
	for pos = 0 : L-1
		tot%=3
		if evenxs
			z3Flag[ind] |= (tot << (2*pos))
		else
			# edge to the left is invertible
			z3Flag[ind] |= (3 << (2*pos))
		end
		a=(state >> (2*pos)) & 3
		if a==0
			evenxs = ! evenxs
		elseif a==2
			tot+=1
		elseif a==3
			tot+=2
		end
	end
end

function stringFromZ3(z3Labels,L)
	s=""
	for i in 1 : L
		s*=revZ3Mapping[(z3Labels&3)]
		z3Labels>>=2
	end
	return s
end

println("preparing...")
for i = 1 : 4^L
	setZ3Flag!(z3Flag_,i,flag_,L)
end

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


#=
Zipper stuff
=#

longFlag_ = zeros(Int32,4^(L+2))
longZ3Flag_ = zeros(Int64,4^(L+2))
longEdgeState_ = zeros(Int64,4^(L+2))
longRevEdgeState_ = zeros(Int64,8^(L+2))

# println("preparing...")
# for i = 1 : 4^(L+2)
# 	setFlag!(longFlag_,i,L+2)
# end
# for i = 1 : 4^(L+2)
# 	setZ3Flag!(longZ3Flag_,i,longFlag_,L+2)
# end
# for i = 1 : 4^(L+2)
# 	setEdgeState!(longEdgeState_,longRevEdgeState_,i,longFlag_,L+2)
# end


# for i = 1 : 4^(L+2)
# 	if mainFlag(longFlag_,i,L+2) == 0
# 		println("state: ", stringFromIndex(i,L+2))
# 		println("edge:  ", stringFromEdgeState(longEdgeState_[i],L+2))
# 		# println("z3:    ", stringFromZ3(z3Flag_[i],L))
# 		println()
# 	end
# end


const ζinvLabels = [(0, 0, 0, 0), (0, 0, 1, 0)]

const ξLabels = [(0, 0, 0, 1), (0, 0, 0, 2), (0, 0, 0, 3), (0, 0, 1, 1), (0, 0, 1, 2),
(0, 0, 1, 3), (0, 1, 1, 0), (0, 2, 3, 0), (0, 3, 2, 0)]

const xLabels = [(0, 1, 1, 1), (0, 2, 3, 3), (0, 3, 2, 2), (1, 1, 2, 2), (1, 2, 1, 1),
(1, 3, 3, 3), (2, 1, 3, 3), (2, 2, 2, 2), (2, 3, 1, 1)]

const zLabels = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 2), (0, 2, 2, 1), (0, 3, 1, 3),
(0, 3, 3, 1), (1, 1, 1, 3), (1, 1, 3, 1), (1, 2, 2, 3), (1, 2, 3, 2),
(1, 3, 1, 2), (1, 3, 2, 1), (2, 1, 1, 2), (2, 1, 2, 1), (2, 2, 1, 3),
(2, 2, 3, 1), (2, 3, 2, 3), (2, 3, 3, 2)]

const y1Labels = [(0, 1, 1, 2), (0, 1, 2, 1), (0, 1, 3, 3), (0, 2, 1, 3), (0, 2, 2, 2),
(0, 2, 3, 1), (0, 3, 1, 1), (0, 3, 2, 3), (0, 3, 3, 2), (1, 1, 1, 1),
(1, 1, 2, 3), (1, 1, 3, 2), (1, 2, 1, 2), (1, 2, 2, 1), (1, 2, 3, 3),
(1, 3, 1, 3), (1, 3, 2, 2), (1, 3, 3, 1), (2, 1, 1, 3), (2, 1, 2, 2),
(2, 1, 3, 1), (2, 2, 1, 1), (2, 2, 2, 3), (2, 2, 3, 2), (2, 3, 1, 2),
(2, 3, 2, 1), (2, 3, 3, 3)]

const y2Labels = [(0, 1, 1, 3), (0, 1, 2, 2), (0, 1, 3, 1), (0, 2, 1, 1), (0, 2, 2, 3),
(0, 2, 3, 2), (0, 3, 1, 2), (0, 3, 2, 1), (0, 3, 3, 3), (1, 1, 1, 2),
(1, 1, 2, 1), (1, 1, 3, 3), (1, 2, 1, 3), (1, 2, 2, 2), (1, 2, 3, 1),
(1, 3, 1, 1), (1, 3, 2, 3), (1, 3, 3, 2), (2, 1, 1, 1), (2, 1, 2, 3),
(2, 1, 3, 2), (2, 2, 1, 2), (2, 2, 2, 1), (2, 2, 3, 3), (2, 3, 1, 3),
(2, 3, 2, 2), (2, 3, 3, 1)]

function Inv!(s)
	if s==2
		return 3
	elseif s==3
		return 2
	else
		return s
	end
end

# Twisted Z3 charge
function Q!(s)
	if s==0
		return 0
	else
		return s-1
	end
end

function Attach!(C,B)
	for ind = 1 : 4^(L+2)
		C[ind] = 0
	end
	for ind = 1 : 4^L
		if mainFlag(flag_,ind,L) != 0
			continue
		end
		state = stateFromInd(ind) << 2
		if (isodd(trailingXs(state))) # start label is 1
			ni = newInd(state,L+2,sXX,L+2)
			# if mainFlag(longFlag_,ni,L+2) == 0
				C[ni] += B[ind]
			# end
		else
			ni = newInd(state,L+2,sXX,L+2)
			# if mainFlag(longFlag_,ni,L+2) == 0
				C[ni] += 1/ζ * B[ind]
			# end
			ni = newInd(state,L+2,s00,L+2)
			# if mainFlag(longFlag_,ni,L+2) == 0
				C[ni] += ξ * B[ind]
			# end
			ni = newInd(state,L+2,sPM,L+2)
			# if mainFlag(longFlag_,ni,L+2) == 0
				C[ni] += ξ * B[ind]
			# end
			ni = newInd(state,L+2,sMP,L+2)
			# if mainFlag(longFlag_,ni,L+2) == 0
				C[ni] += ξ * B[ind]
			# end
			# C[newInd(state,L+2,sXX,L+2)] += 1/ζ * B[ind]
			# C[newInd(state,L+2,s00,L+2)] += ξ * B[ind]
			# C[newInd(state,L+2,sPM,L+2)] += ξ * B[ind]
			# C[newInd(state,L+2,sMP,L+2)] += ξ * B[ind]
		end
	end
end

function ZipF!(z3, s1, s2, s3, s4)
	if z3==3
		# first edge is invertible
		if s1==s3==x && s2==s4
			return 1
		else
			return 0
		end
	end

	if ((Q!(s1)+Q!(s2)-Q!(s3)-Q!(s4))%3) != 0 || ((s1!=0)+(s2!=0)+(s3!=0)+(s4!=0)) != 0
		# s1,s2 and s3,s4 must get to the same third edge
		return 0
	# elseif s1==s2==0
	# 	if s3!=Inv!(s4)
	# 		return 0
	# 	end
	# elseif (s1==0 || s2==0) && (s1!=s3 || s2!=s4)
	# 	return 0
	# elseif s3==s4==0
	# 	if s1!=Inv!(s2)
	# 		return 0
	# 	end
	# elseif (s3==0 || s4==0) && (s1!=s3 || s2!=s4)
	# 	return 0
	end

	label = (z3, s1, s2, s3)
	if label in ζinvLabels
		return 1/ζ
	elseif label in ξLabels
		return ξ
	elseif label in xLabels
		return x
	elseif label in zLabels
		return z
	elseif label in y1Labels
		return y1
	elseif label in y2Labels
		return y2
	else
		return 1
	end
end

function Zip!(C,B,i,longZ3Flag)
	for ind = 1 : 4^(L+2)
		C[ind] = 0
	end
	for ind = 1 : 4^(L+2)
		# if mainFlag(longFlag_,ind,L+2) != 0
		# 	continue
		# end
		z3 = (longZ3Flag[ind] >> 2*(i-1)) & 3
		state = stateFromInd(ind,L+2)
		# if (ind == 4^L)
		# 	println(i)
		# 	println(stringFromIndex(ind,L+2))
		# 	println(localStatePair(state,3,L+2))
		# 	println()
		# end
		s1,s2 = localStatePair(state,i,L+2)
		for s3 = 0 : 3
			for s4 = 0 : 3
				ni = newInd(state,i,(s3,s4),L+2)
				if mainFlag(longFlag_,ind,L+2) == 0
					C[ni] += ZipF!(z3,s1,s2,s3,s4) * B[ind]
				# if (ind == 4^L)
					edgeState = longEdgeState_[ind]
					j = i
					e1 = (edgeState>>(3*(j-1)))&7
					j = nextSite(j,L+2)
					e2 = (edgeState>>(3*(j-1)))&7
					j = nextSite(j,L+2)
					e3 = (edgeState>>(3*(j-1)))&7
					if e1<4
						e4=e1+3
					else
						if s3==0
							e4=e1-3
						else
							e4=4+((e1+s3-2)%3)
						end
					end
					if (ZipF!(z3,s1,s2,s3,s4) != 0) && (FSymbol!(4,e1,4,e3,e2,e4) != ZipF!(z3,s1,s2,s3,s4))
						# println(i)
						println(stringFromIndex(ind,L+2)," at position ",i)
						# println(revmapping[s4])
						println("e1...e4: ",edgeRevmapping[e1],edgeRevmapping[e2],edgeRevmapping[e3],edgeRevmapping[e4])
						# println(4,e1,4,e3,e2,e4)
						println("FSymbol(",4,e1,4,e3,e2,e4,") = ", FSymbol!(4,e1,4,e3,e2,e4))
						# println(edgeState[nextSite(i,L+2)])
						println(stringFromEdgeState(edgeState,L+2))
						# println(z3,s1,s2,s3,s4)
						println("ZipF(",z3,s1,s2,s3,s4,") = ",ZipF!(z3,s1,s2,s3,s4))
						println()
					end
				# end
				end

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
	for ind = 1 : 4^(L+2)
		if mainFlag(longFlag_,ind,L+2) != 0
			continue
		end
		state = stateFromInd(ind,L+2)
		ni = state&(4^L-1)
		if ni==0
			ni=4^L
		end
		if mainFlag(flag_,ni,L) != 0
			continue
		end
		# println()
		# println(state)
		# println(ni)
		# println()

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

# attach = LinearMap((C,B)->Attach!(C,B),4^(L+2),4^L,ismutating=true,issymmetric=false,isposdef=false)
# ρ = attach
#
# for i = 1 : L
# 	zip = LinearMap((C,B)->Zip!(C,B,i,longZ3Flag_),4^(L+2),ismutating=true,issymmetric=false,isposdef=false)
# 	global ρ = zip * ρ
# end

# t = LinearMap((C,B)->Tfunc!(C,B,L+2,false),4^(L+2),ismutating=true,issymmetric=false,isposdef=false)
# ρ = t * ρ

# detach = LinearMap((C,B)->Detach!(C,B),4^L,4^(L+2),ismutating=true,issymmetric=false,isposdef=false)
# ρ = detach * ρ

# println(Matrix(adjoint(v) * ρ * v))

# println(norm(Matrix(ρ)))
# e,v = eigen(Matrix(ρ))
# println(e)
# println(Matrix(ρ))


# edgeState = longEdgeState_
# flag = longFlag_


# str = "000"
# ind = indexFromString(str,inL)

# for ind = 1 : 4^(L+2)
#
# 	inL = L+2
# 	edgeState = longEdgeState_
# 	flag = longFlag_
#
# 	# if mainFlag(flag,ind,inL)==0
#
# 		str = stringFromIndex(ind,inL)
# 		println()
# 		println(str)
#
# 		println(str)
# 		# , " = ", stringFromEdgeState(edgeState[ind],inL))
# 		testV = zeros(4^inL)
# 		testV[ind] = 1
# 		println("mainFlag: ", mainFlag(flag,ind,L)==0)
#
# 		global zip = LinearMap((C,B)->Zip!(C,B,1,longZ3Flag_),4^(L+2),ismutating=true,issymmetric=false,isposdef=false)
# 		testU = zip * testV
#
# 		# testU = attach * testV
#
# 		outL = L+2
# 		flag = longFlag_
# 		# edgeState = edgeState_
# 		edgeState = longEdgeState_
# 		for ind = 1 : 4^outL
# 			if testU[ind] != 0
# 				# if mainFlag(flag,ind,outL) == 0
# 					println(stringFromIndex(ind,outL), " has value ", testU[ind])
# 					# println(stringFromIndex(ind,outL), " = ", stringFromEdgeState(edgeState[ind],outL), " has value ", testU[ind])
# 				# else
# 				# 	println("bad flag: ", ind, " = ", stringFromIndex(ind,outL))
# 				# end
# 			end
# 		end
#
# 	# end
# end



# println()
# smallH = Matrix(diagm(e))
# smallT = Matrix(adjoint(v)*T*v)
# smallρ = Matrix(adjoint(v)*ρ*v)
# println(smallρ)
# smalle,smallv = eigen(smallH+smallT+smallρ)
# Hs = real(diag(adjoint(smallv)*smallH*smallv))
# Ts = diag(adjoint(smallv)*smallT*smallv)
# Ps = real(map(log,Ts)/(2*π*im))*L
# ρs = diag(adjoint(smallv)*smallρ*smallv)
# HPρs = hcat(Hs,Ps,ρs)
# HPρs = sort([HPρs[i,:] for i in 1:size(HPρs, 1)])
# s=""
# s*=string(HPρs[1][1])
#
# print(MathematicaMatrix(HPρs))


# function Zip!(C,B,flag,z3Flag)
# 	for ind = 1 : 4^(L+2)
# 		state = stateFromIndLong(ind)
# 		for i = 1 : L
# 			sp=localStatePair(state,i)
# 			if sp==sXX  && isρ1ρ(flag,ind,i)
# 				C[newIndLong(state,i,sPM)] -= ξ * y1 * B[ind]
# 				C[newIndLong(state,i,sMP)] -= ξ * y2 * B[ind]
# 				C[newIndLong(state,i,s00)] -= ξ * x * B[ind]
# 			elseif sp==sPM
# 				C[newIndLong(state,i,sXX)] -= y1 * ξ * B[ind]
# 				C[newIndLong(state,i,sMP)] -= y1 * y2 * B[ind]
# 				C[newIndLong(state,i,s00)] -= y1 * x * B[ind]
# 			elseif sp==sMP
# 				C[newIndLong(state,i,sXX)] -= y2 * ξ * B[ind]
# 				C[newIndLong(state,i,sPM)] -= y2 * y1 * B[ind]
# 				C[newIndLong(state,i,s00)] -= y2 * x * B[ind]
# 			elseif sp==s00
# 				C[newIndLong(state,i,sXX)] -= x * ξ * B[ind]
# 				C[newIndLong(state,i,sPM)] -= x * y1 * B[ind]
# 				C[newIndLong(state,i,sMP)] -= x * y2 * B[ind]
# 			elseif sp==s0P
# 				C[newIndLong(state,i,sP0)] -= y1 * y2 * B[ind]
# 				C[newIndLong(state,i,sMM)] -= y1 * z * B[ind]
# 			elseif sp==sP0
# 				C[newIndLong(state,i,s0P)] -= y2 * y1 * B[ind]
# 				C[newIndLong(state,i,sMM)] -= y2 * z * B[ind]
# 			elseif sp==sMM
# 				C[newIndLong(state,i,s0P)] -= z * y1 * B[ind]
# 				C[newIndLong(state,i,sP0)] -= z * y2 * B[ind]
# 			elseif sp==s0M
# 				C[newIndLong(state,i,sM0)] -= y2 * y1 * B[ind]
# 				C[newIndLong(state,i,sPP)] -= y2 * z * B[ind]
# 			elseif sp==sM0
# 				C[newIndLong(state,i,s0M)] -= y1 * y2 * B[ind]
# 				C[newIndLong(state,i,sPP)] -= y1 * z * B[ind]
# 			elseif sp==sPP
# 				C[newIndLong(state,i,s0M)] -= z * y2 * B[ind]
# 				C[newIndLong(state,i,sM0)] -= z * y1 * B[ind]
# 			else
# 				return
# 			end
# 		end
# 	end
# end
