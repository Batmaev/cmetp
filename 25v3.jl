### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 59db45bf-c9f7-4989-b92a-f0862bd7a474
begin
	
	import Pkg
	Pkg.activate(mktempdir())
	
	Pkg.add("PlutoUI")
	using PlutoUI
	
	Pkg.add("Plots")
	using Plots
	plotly()
	
end

# ╔═╡ efbb1b00-1e4b-11ec-32e0-3b41bf2b8b36
md"""

### Задание 2, № 5 по книге
###### Постановка задачи:

Уравнение:

``
\displaystyle
\frac{\partial u}{\partial t} + \frac{\partial F}{\partial x}
= \nu \frac{\partial^2 u}{\partial x^2},
\qquad F = \frac{u^2}{2}
``

Начальное условие: 

a) ``\; \displaystyle 
u|_{t=0} = \operatorname{HeavisideTheta}(-x) 
``

b) ``\; \displaystyle
u|_{t=0} = \exp\!\left(
	-\frac{(x-x_0)^2}{d^2}
\right)
``

Краевые условия:

`` u(t, a) \equiv u(0, a) \qquad u(t,b) \equiv u(0,b) ``

Числа Куранта:

`` \displaystyle \mu = \frac{\nu \tau}{h^2}, 
\qquad r = \frac{\tau}{h}``


"""

# ╔═╡ 0656c12b-d327-47dc-9107-2ed74f1f2623
md"""
###### Зададим параметры и начальные функции:
"""

# ╔═╡ 9a2596d1-493a-45f0-b8bf-7f379d667d58
begin
	a = -1           # Границы промежутка по пространству
	b = 1
	
	T = 1            # Граница промежутка по времени

	Nx = 1000        # Точки нумеруются от 1 до nX+1
	Nt = 5000
	
	τ = T / Nt
	
	h = (b - a) / Nx
	
	
	ν = 0.01 		 # Коэффициент в уравнении
	
	r = τ / h        # Гиперболическое число Куранта
	
	μ = ν * τ / h^2  # Параболическое число Куранта
end

# ╔═╡ d7f30bc9-a783-4b0e-ad1e-d60b166f3f25
u0a(x) = x > 0 ? 0 : 1

# ╔═╡ 2ad832ed-d147-4fc9-a404-a0c4b4e13871
begin
	x0 = 0
	d = 0.1
	
	u0b(x) = exp(-(x-x0)^2 / d^2)
end

# ╔═╡ cf19dfb5-2180-4a07-87e2-113d58ff81ab
# Выбираем начальную функцию из 2 вариантов
u0 = u0b

# ╔═╡ 8730a761-5ed2-49ef-85e3-738d9bdbff40
md"""

###### Используем схему Мак-Кормака:

``\displaystyle
r = \frac{\tau}{h},
\qquad
\mu = \frac{\nu \tau}{h^2},
\qquad
F = \frac{u^2}{2}
``

1). Промежуточный временной слой:

``\qquad 
\bar u_j = u_j - r(F_{j+1} - F_j) + \mu (u_{j+1} - 2 u_j + u_{j-1})
``

2). Следующий временной слой:

``\qquad
\hat u_j = \frac{1}{2}(u_j + \bar u_j) - \frac{r}{2}(\bar F_j - \bar F_{j-1}) + \frac{\mu}{2} (\bar u_{j+1} - 2 \bar u_j + \bar u_{j+1}) ``

"""

# ╔═╡ 974b08df-41e0-4917-af66-095c1c6a4032
md"""
###### Устойчивость:

Чтобы схема Мак-Кормака была устойчивой, необходимо одновременное выполненение двух условий:

``\displaystyle
\mu \leqslant \frac{1}{2}
\qquad\quad
r \cdot \operatorname{max} |u| \leqslant 1
``
"""

# ╔═╡ bb90c7b3-3bbd-42d3-bc05-4545c5cfa8dc
md"""
У нас:

``\mu = `` $μ

``\operatorname{max}|u| \overset{?}{=} 1 ``

``r = `` $r

"""

# ╔═╡ 204c1031-52e5-4a7d-97e1-9e2054c425dd
md"""
###### Точки возле границ:

Если обозначить настоящие узлы как `*`, а промежуточные -- как `⋅`, то шаблон можно изобразить так:

```
        *
        |
    ⋅ - ⋅ - ⋅
    |   |   |
* - * - * - * - * 
```

Кажется, что схема Мак-Кормака не позволяет определить первые две и последние две точки. Однако мы можем при помощи краевых условий доопределить полностью значения сначала на промежуточном временном слое, а затем и на настоящем временном слое.  
"""

# ╔═╡ 64311ba3-9ce5-49ca-8667-2b40459ad384
md"""
###### Код:
"""

# ╔═╡ 73e25696-0c66-4b97-b156-9a6b5e93338f
# Инициализация
begin
	xlist = a:h:b
	
	# Здесь будет храниться решение
	u  = Array{Float64}(undef, Nt+1, Nx+1)
	
	# Первый временной слой
	u[1, :] = u0.(xlist)
	
	# Промежуточный временной слой (будет перезаписываться на каждой итерации)
	uBar = Array{Float64}(undef, Nx+1)
	
	# F = U^2/2, где в качестве U может быть 
	# как настоящий, так и промежуточный временной слой
	F = [U^2/2 for U ∈ u[1, :]]
end;

# ╔═╡ c63ad377-2157-4ea8-8121-e4414bff4a6a
# Основной цикл

for n = 1 : Nt
	
	F[:] = [U^2/2 for U ∈ u[n, :]]
	
	# Промежуточный временной слой
	for j = 2 : Nx
		uBar[j] = u[n,j] - r*(F[j+1] - F[j]) + μ*(u[n,j+1] - 2u[n,j] + u[n,j-1])
	end
	uBar[1] = u0(a)
	uBar[Nx+1] = u0(b)
	
	F[:] = [U^2/2 for U ∈ uBar]
	
	
	# Новый временной слой
	for j = 2 : Nx
		u[n+1,j] = 1/2*(u[n,j] + uBar[j]) 
			     - r/2*(F[j] - F[j-1]) 
			     + μ/2 * (uBar[j+1] - 2uBar[j] + uBar[j-1])
	end
	u[n+1, 1] = u0(a)
	u[n+1, Nx+1] = u0(b)
	
end

# ╔═╡ e35eb3f0-6aa8-41dc-89db-3710df7e56ec
@bind t Slider(1 : Nt)

# ╔═╡ 44e78562-f35d-41e3-9fa9-b9fd0a16b800
plot(
	xlist, u[t, :],
	ylims=(0, 1.1),
	title = "t = $(round((t*τ), digits=3))",
	label="",
)

# ╔═╡ Cell order:
# ╠═59db45bf-c9f7-4989-b92a-f0862bd7a474
# ╟─efbb1b00-1e4b-11ec-32e0-3b41bf2b8b36
# ╟─0656c12b-d327-47dc-9107-2ed74f1f2623
# ╠═9a2596d1-493a-45f0-b8bf-7f379d667d58
# ╠═d7f30bc9-a783-4b0e-ad1e-d60b166f3f25
# ╠═2ad832ed-d147-4fc9-a404-a0c4b4e13871
# ╠═cf19dfb5-2180-4a07-87e2-113d58ff81ab
# ╟─8730a761-5ed2-49ef-85e3-738d9bdbff40
# ╟─974b08df-41e0-4917-af66-095c1c6a4032
# ╟─bb90c7b3-3bbd-42d3-bc05-4545c5cfa8dc
# ╟─204c1031-52e5-4a7d-97e1-9e2054c425dd
# ╟─64311ba3-9ce5-49ca-8667-2b40459ad384
# ╠═73e25696-0c66-4b97-b156-9a6b5e93338f
# ╠═c63ad377-2157-4ea8-8121-e4414bff4a6a
# ╟─e35eb3f0-6aa8-41dc-89db-3710df7e56ec
# ╠═44e78562-f35d-41e3-9fa9-b9fd0a16b800
