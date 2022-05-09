### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 63def3e5-bb09-4b9f-8eff-5f71bceaf2cd
using PlutoUI

# ╔═╡ fc59a3b3-4343-46cf-bbec-858e049aa32b
import LinearAlgebra: norm

# ╔═╡ abb35e43-c186-4066-9e74-fd06b4b35f45
TableOfContents()

# ╔═╡ 2977679a-a1c4-48ac-a009-7f74516babce
md"# 🧾 Background

This `Pluto.jl` notebook is a simple runnable demonstration of [space-filling curves](https://en.wikipedia.org/wiki/Space-filling_curve).

I recommend the following book

```
Ventrella, Jeffrey. Brainfilling curves-a fractal bestiary. Lulu. com, 2012.
```

as a good starting point for further reading. A great deal of inspiration for this notebook was drawn from that excellent freely avaible resource. You can [read it online here](http://www.brainfillingcurves.com/)

The goal of this notebook is to show some runnable code as well as to show off `Pluto.jl` notebooks which are always enjoyable to work with.
"

# ╔═╡ fcdb08f3-f219-423d-912e-4085eeff45da
md"## Preparing our tools

Never work without the right equipment. Let's define some structs that will make our life easier, and make the code clearer."

# ╔═╡ 17ab22e0-add4-455b-8896-0184f6fe2eec
mutable struct Line
	rotation::Float64
	length::Float64
	mirror::Vector{Bool}
end

# ╔═╡ 9aad0ae4-8901-42bf-bb1c-8e92129b64c8
md"Let's implement copy for this new struct `line` so that we can perform copy operation"

# ╔═╡ 11d2de71-790a-4397-a0e1-ad1ebc40b82b
Base.copy(s::Line) = Line(s.rotation, s.length, copy(s.mirror))

# ╔═╡ 97b235f9-2b55-451d-b6c0-e40de4e91c4d
md"To make our code a bit nicer to write, lets write a simple macro to create line objects. Instead of `Line(30, 2, [1,0])` one can write `@Line 30 2 twist`. For more details on this notation and interpretation, see [the following section of the online resource](https://archive.org/details/BrainfillingCurves-AFractalBestiary/page/n47/mode/2up?view=theater)"

# ╔═╡ 33890185-0f3c-4460-a37c-41e5ffbf5c1f
macro Line(angle, length, args...)
	mirror = [false,false]
	for arg in args
		if arg==:twist mirror[2]=true  
		elseif arg==:reverse mirror[1]=true end		
	end
	return Line(eval(angle), eval(length), mirror)
end;

# ╔═╡ 0f80fc16-e457-4d48-8116-ab8e854a6183
md"Let's also create a struct for a 'template' that will be used for each iteration. More or less a glorified list of `line` structs"

# ╔═╡ c86b7a75-4448-4313-a962-037bef8bc033
md"We need a simple utility function that helps to determine where a series of lines will end up."

# ╔═╡ 42321e29-a4f6-4a8a-a25e-ec445176c6ff
function computeEndPoint(lines::Vector{Line})::Vector{Float64}
	return sum(lines) do line
		return [cos(deg2rad(line.rotation)), sin(deg2rad(line.rotation))]*line.length 
	end
end

# ╔═╡ ca1776ad-3bc0-4a03-ac74-526e1f1e0950
struct LineTemplate
	templateLines::Vector{Line}
	endpoint::Vector{Float64}
	span::Float64
	rotation::Float64
	LineTemplate(template) = new(template, computeEndPoint(template), norm(computeEndPoint(template)), rad2deg(atan(computeEndPoint(template)[2],computeEndPoint(template)[1])))
end

# ╔═╡ 89d8e415-825f-4db8-a493-155c7304d5c3
function computeBBox(lines::Vector{Line})::NTuple{4, Float64}
	minX, minY, maxX, maxY, x, y = (0., 0., 0., 0., 0., 0.)
	for line in lines
		x += cos(deg2rad(line.rotation))*line.length
		if (x > maxX) maxX = x elseif (x < minX) minX = x end		

		y += sin(deg2rad(line.rotation))*line.length
		if (y > maxY) maxY = y elseif (y < minY) minY = y end
	end	
	# return minX, minY, width, height
	return minX, minY, maxX-minX, maxY-minY 
end

# ╔═╡ bb9b76bc-9e79-4fb3-bf98-592c14a38335
md"## Defining iteration functions

Now we need to define the core of the process: the iteration utilities. 

Let's start with an interface function. It requires a `seed` (a starting shape to start iterating on with the  `template`) and a `template` which is what each `line` will be replaced with recursively.

> Note that in the most basic case you choose the template itself to be the seed.

And of course a number of iterations must be provided.
"

# ╔═╡ fb205fe8-ad82-4c9a-8312-8d8b8554ea46
function iter_curve(out::Vector{Line}, sourceline::Line, template::LineTemplate, depth::Int64)	
	# if the sourceline is reversed, then we must go through template in reverse
	templatelist =(sourceline.mirror[1] ? reverse(template.templateLines) : template.templateLines)
	for templateLine in templatelist	
		line = copy(templateLine)  # dont want to affect original template
		scalefac = line.length/template.span		
		line.length = scalefac * sourceline.length
		
		line.rotation -= template.rotation
		if sourceline.mirror[2]
			line.mirror[2] = !line.mirror[2]
			line.rotation *= -1
		end
		if sourceline.mirror[1]
			line.mirror[1] = !line.mirror[1]  
			line.rotation *= -1
		end
		line.rotation += sourceline.rotation
		
		if depth > 0			
			iter_curve(out, line, template, depth-1)
		else
			push!(out, line)
		end
	end	
end

# ╔═╡ a053a0f4-177e-4e82-b804-69d71fa7ba98
function computeSFC(seed::Vector{Line}, template::LineTemplate, iterations::Int64)::Vector{Line}
	@assert iterations >= 0 "Number of iterations must be positive"

	curve::Vector{Line} = []
	
	for line in seed
		if iterations > 0 
			iter_curve(curve, line, template, iterations-1)
		else
			push!(curve, line)
		end
	end	

	return curve
end

# ╔═╡ 14f161e7-c0b4-48d2-a37c-85ffd5c35dd7
md"## Visualisation

We need to be able to view our curves. `Pluto.jl` will show the text/html representation of the object, which we can define ourselves! So we draw our curves directly on a SVG canvas.

Below are the main colours used for the plot, which you may want to tune for your curve figures."

# ╔═╡ f166b2db-3279-4d63-98fd-cc87cf02c0d5
begin
	linecolor = "5CDB95"
	pivotcolor = "05386B"
	backgroundcolor = "EDF5E1"
end;

# ╔═╡ ee0bc1fe-ca91-4409-ac1d-bdda754f4a65
md"We can define any group of valid SVG elements to make up our 'line' used for visualisation. Note that it would be good to have it go from `0 < x < 1` and have it look somewhat line-like"

# ╔═╡ 2bae6d42-7506-43ca-9947-304a5436efe0
lineshape = """
<line x1="0.08" y1="0" x2="0.92" y2="0" stroke="#$(linecolor)" stroke-width="0.1" stroke-linecap="round"/>
<line x1="0.8" y1="0.2" x2="0.92" y2="0" stroke="#$(linecolor)" stroke-width="0.1" stroke-linecap="round"/>""";

# ╔═╡ e0429a1b-b96c-4f03-8e69-1efc2c87054a
md">But you can also go crazy any make a line out of e.g. emojis 😜. Replace `lineshape` with `lineshapeEmoji` in the `Base.show` function below to try it."

# ╔═╡ bc8f2c5a-ba8d-42a3-8d43-21c6a3b6a0a9
lineshapeEmoji = """
<text text-anchor="middle" alignment-baseline="middle"font-size="0.03em" transform="translate(0.45,0.05) rotate(45)">💉</text>""";

# ╔═╡ c62a37c7-9b4d-4c04-bc58-140ad69c94dd
md"The method below is not elegant code, but it gets the job done. You can dig around a bit in the method below if you want more fine grained control over the plots (change shapes, sizes etc)

We add a simple ugly button that appears on hover that will download the SVG file. Then you can open this in your browser to get a fullscreen view or import it to a vector graphics editor.."

# ╔═╡ acaf5695-f92b-485f-8ed8-461adab36c3e
begin
	function Base.show(io::IO, ::MIME"text/html", obj::Vector{Line})
		(minX, minY, width, height) = computeBBox(obj)
		svg = """
		<span class="SFC-fig" style="position:relative">
		<button onclick="const svg = this.nextElementSibling.outerHTML; const blob = new Blob([svg.toString()]); const element = document.createElement('a'); element.download = 'space-filling-curve.svg'; element.href = window.URL.createObjectURL(blob); element.click(); element.remove()"	
		type="button"> Download SVG </button>		
		<svg viewBox="$(minX-0.3) $(minY-0.3) $(width+0.6) $(height+0.6)" style="background:#$(backgroundcolor);transform:scaleY(-1)" xmlns="http://www.w3.org/2000/svg"> 
		<defs>
			<g id="curve-line"> $(lineshape) </g>
		</defs>
		<circle cx="0" cy="0" r="0.03" fill="black" opacity="0.3"/>
		"""
		lastcoord = [0,0]		
		for line in obj
			# add to svg
			scalefacX = line.mirror[1] ? -line.length : line.length
			scalefacY = line.mirror[2] ? -line.length : line.length
			svg *= """<use href="#curve-line" transform="translate($(lastcoord[1]),$(lastcoord[2])) rotate($(line.rotation)) $(line.mirror[1] ? "translate(1,0)" : "") scale($(scalefacX),$(scalefacY)) "/>
			<circle cx="$(lastcoord[1])" cy="$(lastcoord[2])" r="$(line.length/40)" fill="#$(pivotcolor)"/>"""
			
			lastcoord= lastcoord + [cos(deg2rad(line.rotation)), sin(deg2rad(line.rotation))]*line.length
		end	
		svg *= """</svg>		
		<style>
		.SFC-fig button { display:none; position:absolute; z-index:1 }
		.SFC-fig:hover button { display:block }
		</style>
		</div>"""
	    write(io, svg)
	end
end;

# ╔═╡ 391a4b81-4ca7-4df4-add6-d91094224cc8
md"We ensure that our `LineTemplate` struct gets displayed similar to Vector{Line}"

# ╔═╡ 3dadeb2e-16ca-43b3-a3d9-035f58911d4c
Base.show(io::IO, m::MIME"text/html", obj::LineTemplate) = Base.show(io, m,obj.templateLines)

# ╔═╡ f4be52c2-e6a1-4c5e-89e2-83361abd5262
md"# 🚀 Let's draw some curves!"

# ╔═╡ ad207712-1781-45d8-b462-3025a6c4baa8
md"## Simplest Curves

Let's start with the simplest template we can think of for a curve, two `line` segments connected at 90 degrees. 

Changing the mirror properties of the segments will change the resulting pattern rather dramatically."

# ╔═╡ 11bf8396-aeb5-4dec-a4e2-ff26a2912297
rightangle = LineTemplate([
		@Line 0 1 twist
		@Line 90 1 
])

# ╔═╡ 25ab6b0b-33fe-437c-bcf0-cdbac6816a9a
md"Try playing with the `twist` and `reverse` values in the template to see how the curve changes!"

# ╔═╡ 5a045687-2a71-44d0-995b-222950fd203d
@bind rightangleorder Slider(0:10, default=9, show_value=true)

# ╔═╡ e655ca64-f509-459c-b8bf-31811a622425
computeSFC(rightangle.templateLines, rightangle, rightangleorder)

# ╔═╡ 3a47428b-d5f7-4540-83f8-139c6b7359f3
md"## Koch Curve

Let's continue with a classic: the Koch Curve"

# ╔═╡ 0beafe56-a4c6-4ade-8954-43998ba5c2a0
kochcurve = LineTemplate([
		@Line 0 1
		@Line 60 1
		@Line -60 1
		@Line 0 1
])

# ╔═╡ dcd0f0c7-9a32-46c7-b932-567dcbde6e3e
@bind kochorder Slider(0:5, default=4, show_value=true)

# ╔═╡ 279b3826-c78b-4497-945b-5e3f767f854c
computeSFC(kochcurve.templateLines, kochcurve, kochorder)

# ╔═╡ 90375842-b83a-4062-8fee-7ba768d6ca9e
md"We are not required to use the template as the base curve at all, so we can try some interesting seed geometries! Lets try the koch curve template starting from a square (check seed by setting iterations/order to 0)"

# ╔═╡ 2bff28ae-87e8-4bd9-9fa6-e6e143872c3f
computeSFC([
		@Line 0 1
		@Line 90 1
		@Line 180 1
		@Line 270 1 
], kochcurve, kochorder)

# ╔═╡ 9c37f943-e7a2-4ad8-8916-c00938ffdbed
md"## Holiday Tree by Jeffrey Ventrella

This is a pretty nice curve, and requires the use of orientations of our `line` strucs."

# ╔═╡ 85cad6a9-2b14-4bed-9c15-81145a96c3ed
holidayTree = LineTemplate([
		@Line 30 sqrt(3) twist 
		@Line 120 1  
		@Line 0 1 
		@Line -120 1 
		@Line -30 sqrt(3) twist 
])

# ╔═╡ 796a211a-bb0d-4f54-a1f1-18d94ffa9439
@bind holidayorder Slider(0:5, default=4, show_value=true)

# ╔═╡ 77ba6866-7a8c-4451-8705-c9ccbddecccd
computeSFC(holidayTree.templateLines, holidayTree, holidayorder)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═63def3e5-bb09-4b9f-8eff-5f71bceaf2cd
# ╠═fc59a3b3-4343-46cf-bbec-858e049aa32b
# ╠═abb35e43-c186-4066-9e74-fd06b4b35f45
# ╟─2977679a-a1c4-48ac-a009-7f74516babce
# ╟─fcdb08f3-f219-423d-912e-4085eeff45da
# ╠═17ab22e0-add4-455b-8896-0184f6fe2eec
# ╟─9aad0ae4-8901-42bf-bb1c-8e92129b64c8
# ╠═11d2de71-790a-4397-a0e1-ad1ebc40b82b
# ╟─97b235f9-2b55-451d-b6c0-e40de4e91c4d
# ╠═33890185-0f3c-4460-a37c-41e5ffbf5c1f
# ╟─0f80fc16-e457-4d48-8116-ab8e854a6183
# ╠═ca1776ad-3bc0-4a03-ac74-526e1f1e0950
# ╟─c86b7a75-4448-4313-a962-037bef8bc033
# ╠═42321e29-a4f6-4a8a-a25e-ec445176c6ff
# ╠═89d8e415-825f-4db8-a493-155c7304d5c3
# ╟─bb9b76bc-9e79-4fb3-bf98-592c14a38335
# ╠═a053a0f4-177e-4e82-b804-69d71fa7ba98
# ╠═fb205fe8-ad82-4c9a-8312-8d8b8554ea46
# ╟─14f161e7-c0b4-48d2-a37c-85ffd5c35dd7
# ╠═f166b2db-3279-4d63-98fd-cc87cf02c0d5
# ╟─ee0bc1fe-ca91-4409-ac1d-bdda754f4a65
# ╠═2bae6d42-7506-43ca-9947-304a5436efe0
# ╟─e0429a1b-b96c-4f03-8e69-1efc2c87054a
# ╠═bc8f2c5a-ba8d-42a3-8d43-21c6a3b6a0a9
# ╟─c62a37c7-9b4d-4c04-bc58-140ad69c94dd
# ╠═acaf5695-f92b-485f-8ed8-461adab36c3e
# ╟─391a4b81-4ca7-4df4-add6-d91094224cc8
# ╠═3dadeb2e-16ca-43b3-a3d9-035f58911d4c
# ╟─f4be52c2-e6a1-4c5e-89e2-83361abd5262
# ╟─ad207712-1781-45d8-b462-3025a6c4baa8
# ╠═11bf8396-aeb5-4dec-a4e2-ff26a2912297
# ╟─25ab6b0b-33fe-437c-bcf0-cdbac6816a9a
# ╠═5a045687-2a71-44d0-995b-222950fd203d
# ╠═e655ca64-f509-459c-b8bf-31811a622425
# ╟─3a47428b-d5f7-4540-83f8-139c6b7359f3
# ╠═0beafe56-a4c6-4ade-8954-43998ba5c2a0
# ╠═279b3826-c78b-4497-945b-5e3f767f854c
# ╠═dcd0f0c7-9a32-46c7-b932-567dcbde6e3e
# ╟─90375842-b83a-4062-8fee-7ba768d6ca9e
# ╠═2bff28ae-87e8-4bd9-9fa6-e6e143872c3f
# ╟─9c37f943-e7a2-4ad8-8916-c00938ffdbed
# ╠═85cad6a9-2b14-4bed-9c15-81145a96c3ed
# ╠═77ba6866-7a8c-4451-8705-c9ccbddecccd
# ╠═796a211a-bb0d-4f54-a1f1-18d94ffa9439
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
