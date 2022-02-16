### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ 3c9376a6-48cb-41a2-bc8a-2706d5b89a52
begin
	import Pkg
	Pkg.activate(; temp = true)
	
	Pkg.add("Polynomials")
	using Polynomials

	Pkg.add("SpecialPolynomials")
	using SpecialPolynomials

	Pkg.add("FFTW")
	using FFTW

	Pkg.add("Plots")
	using Plots
	plotly()
	theme(:dark)

	Pkg.add("PlutoUI")
	using PlutoUI
end

# ╔═╡ f51b6d70-8dc2-11ec-2f6c-8b52133958a9
md"""
## Задание 6, № 27 по книге
"""

# ╔═╡ 0ca0b86c-c234-482d-bb39-0be8175a421c
md"""
Уравнение Шрёдингера для гармонического осциллятора может быть записано следующим образом:
```math
i \frac{\partial \psi}{\partial t} = 
- \frac12 \frac{\partial^2 \psi}{\partial x^2}
+ \frac{x^2}{2} \psi \tag{1}
```
Требуется найти ``|\psi(t)\rangle``, если:

a) ``\quad |\psi(0)\rangle = |n\rangle`` -- стационарное состояние

b) ``\quad |\psi(0)\rangle = |\alpha\rangle`` -- когерентное состояние
"""

# ╔═╡ 672fece0-3f1b-4f1e-a636-2257b16dd2c8
md"""
###### Сетка
"""

# ╔═╡ db56e6ff-8dd8-456f-8982-7eb11aea7978
md"""
Так как при ``|x| \to \infty`` потенциал неограниченно растёт, волновая функция при больших ``x`` стремится к ``0`` и мы можем ограничиться конечным промежутком ``x \in [-L,\ L]``. 

Введём сетку
```math
x_n = nh - L
,\quad 
h = \frac{2L}{N-1}
,\quad
n = 1 \ldots N-1
\tag{2}
```

"""

# ╔═╡ c2e79179-5bf2-4c79-a293-f3d90d438f95
md"""
#### Метод расщепления по процессам
"""

# ╔═╡ 9eac78d8-cb18-4a32-bd64-830f4d041a75
md"""
Можно разделить шаг по времени продолжительностью ``\tau`` на 2 подшага (тоже продолжительностью ``\tau``). В течение 1 подшага в правой части уравнения (1) будет действовать только первое слагаемое, в течение второго -- только второе.
```math
\begin{align*}
i \frac{\partial \psi}{\partial t} &= - \frac12 \frac{\partial^2 \psi}{\partial x^2} 
\tag{3.1} \\
i \frac{\partial \psi}{\partial t} &= \frac{x^2}{2} \psi \tag{3.2}
\end{align*}
```
Такое приближение на каждом шаге даёт погрешность ``O(\tau^2)``. Более точной оказывается симметричная факторизация, когда второй процесс действует в течение времени ``\tau/2``, затем первый действует в течение времени ``\tau`` и затем снова второй в течение времени ``\tau/2``.

Это с точностью ``O(\tau^3)`` эквивалентно тому, чтобы оба процесса действовали одновременно в течение времени ``\tau``.
"""

# ╔═╡ 82b3aea2-870a-437a-b977-7b94068aaef2
md"""
###### Свободная эволюция
"""

# ╔═╡ 9c20ddb6-19a5-4e90-907f-2c55d1a7a5dc
md"""
Уравнение (2.1) решим с помощью дискретного преобразования Фурье:
```math
\begin{multline}
\psi(t,x_n) = \frac1N \sum_{\nu \in \mathcal N} \widetilde \psi_\nu(t) e^{2\pi i \nu n} 
= \frac1N \sum_{\nu \in \mathcal N} \widetilde \psi_\nu(t) e^{2\pi i \nu (x_n + L) / h}
\end{multline}
```
Здесь ``\nu = \frac1{\lambda}`` -- пространственная частота, ``\mathcal N`` -- множество пространственных частот, а ``n`` связано с ``x_n`` через ``(2)``.

Вычислим производную по ``x_n`` и подставим в уравнение:
```math
i \sum_\nu \dot{\widetilde \psi}_\nu(t) e^{(\dots)} = 
\frac12 \sum_\nu \left( \frac{2\pi \nu }{h} \right)^2 \widetilde \psi_\nu(t) e^{2\pi i \nu (x_n + L) / h}
```
Достаточно обеспечить
```math
i \cdot \dot{\widetilde \psi}_\nu(t) = \frac12 \left( \frac{2\pi \nu }{h} \right)^2 \widetilde \psi_\nu(t)
```
```math
\widetilde \psi_\nu(t + \tau) = \widetilde \psi_\nu(t) \cdot
e^{-\frac{i}{2} \left( \frac{2\pi \nu}{h} \right)^2 \cdot \ \tau}
```

То есть на первом подшаге мы делаем дискретное преобразование Фурье, умножаем коэффициенты на ``e^{-\frac{i}{2} \left( \frac{2\pi \nu}{h} \right)^2 \cdot \ \tau}`` и делаем обратное преобразование Фурье.

"""

# ╔═╡ 139c012b-b1b0-4b93-9739-5ca56230d50e
md"""
###### Множество ``\mathcal N``
Дискретное преобразование Фурье обычно записывают в виде
``\displaystyle
\widetilde \psi_k = \sum_n \psi(x_n) e^{-2\pi i k n / N},
``
то есть в качестве пространственных частот берут ``\nu_k = k/N``. Тогда волновое число ``q = 2\pi \nu \in [0,\ 2\pi)``.

Но физически более правильно брать пространственные частоты такие, что ``q \in [-\pi, \pi)``. Для этого нужно вычесть ``2\pi`` из слишком больших волновых чисел. 

Это не влияет на точки ``x_n``, но влияет на то, как мы интерполируем функцию между ними. Более маленькие частоты вписываются в Теорему Котельникова и не создают высокочастотные осцилляции.

Корректное множество допустимых ``\nu`` можно получить с помощью функции `ftfreq.`
"""

# ╔═╡ 7b45a662-a47d-4b79-8ff9-dcce4254f08e
md"""
###### Эволюция, вызванная потенциалом
"""

# ╔═╡ 0a56188c-5332-4d97-91fe-a69d901d8925
md"""
Уравнение (2.2) решим через оператор эволюции:
```math
\psi(t+\tau, x) = \hat U(\tau) \psi(t, x)
\qquad \hat U(\tau) = e^{-i \hat V \tau}
\qquad \hat V = \frac{x^2}{2}
```
"""

# ╔═╡ e26bff40-4335-43cb-8332-4ef4203821bb
begin
	L  = 10. 		# Будем решать на промежутке [-L, L]
	Nx = 512 		# Число точек по пространств
	h  = 2L / (Nx-1)

	T  = 3π 		# Время моделирования
	τ  = 0.01 		# Шаг по времени
	Nt = ceil(Int, T/τ)
end

# ╔═╡ 18ea867c-43e5-4d7f-af08-0e73662c393b
V(x) = x^2 / 2; 	# Потенциал

# ╔═╡ b7fe7140-8b25-4d7d-86b9-e0436c5af3d9
params = (; L, τ, V);

# ╔═╡ 9f70a119-13de-4739-bef9-cf2e1b6dfefe
function dsolve!(ψ; L, τ, V)
	
	Nx, Nt = size(ψ)
	h = 2L / (Nx - 1)

	U = exp.(-im * τ/2 * V.(-L:h:L))		# Оператор эволюции, вызванной потенциалом

	Ũ = [exp(-im/2 * (2π*ν/h)^2 * τ) for ν ∈ fftfreq(Nx)]
											# Оператор свободной эволюции в Фурье-пространстве

	FT = plan_fft!(ψ[:, 1]) 				# Оператор преобразования Фурье

	for t ∈ 2:Nt
		ψ̃ = @view ψ[:, t]
		ψ̃ .= ψ[:, t-1]
		# Эволюция, вызванная потенциалом, в течение времени τ/2:
		ψ̃ .*= U
		# Свободная эволюция в течение времени τ:
		FT * ψ̃
		ψ̃ .*= Ũ
		FT \ ψ̃
		# Эволюция, вызванная потенциалом, в течение времени τ/2:
		ψ̃ .*= U
	end
end

# ╔═╡ a4a4e044-1956-4efb-8237-b78c67be3174
md"""
###### Стационарное состояние
"""

# ╔═╡ 808fc5bb-8eab-4c2e-9772-ef132761ed17
md"""
```math
|n\rangle = \frac1{\pi^{1/4}\sqrt{n!}\, 2^{n/2}} H_n(x) e^{-x^2/2},
```
где ``H_n(x)`` -- полиномы Эрмита, которые я буду вычислять с помощью пакета SpecialPolynomials.jl.
"""

# ╔═╡ 666e1ffe-335c-4781-8bc0-a695207f03df
Fock(n, x) = 1 / (π^(1/4) * √(factorial(n) * 2^n)) *
		basis(Hermite, n)(x) *
		exp(-x^2 / 2)

# ╔═╡ 97a66676-bd5b-4772-a4ac-3364ea8cb619
ψ₁ = Array{Complex{Float64}}(undef, Nx, Nt);

# ╔═╡ ba6efd83-bd2a-4870-a350-192d7a213748
n = 5;
# Номер стационарного состояния

# ╔═╡ 1ea38beb-bfa5-4f3f-9f58-ce7978f4e087
begin
	ψ₁[:, 1] .= Fock.(n, -L:h:L)
	dsolve!(ψ₁; params...)
end;

# ╔═╡ 8da9463c-a2b6-4e74-82f8-b8c15e2f5bc8
# @bind t₁ Slider(1:Nt, default = 0)
@bind t₁ Clock(0.04, true, false, Nt)

# ╔═╡ aa2528b7-ea47-44ef-9b18-f9a4f93fd046
plot(-L:h:L, 
	[abs.(ψ₁[:, t₁]) angle.(ψ₁[:, t₁])],
	layout = (2,1),
	xlims = (-0.9L, +0.9L),
	ylims = [(-Inf, Inf) (-π, π)],
	labels = ["|ψ|" "arg ψ"],
	ylabel = ["|ψ|" "arg ψ"],
	legend = nothing,
	title = ["|n〉при n = $n и t = $(round((t₁-1)τ, digits = 2))" ""]
)

# ╔═╡ 7566f079-148b-4207-ba77-4f5355b6c394
md"""
При ``t = 0`` волновая функция действительная, положительные числа имеют фазу ``0``, а отрицательные ``\, \pm \pi``.
"""

# ╔═╡ 015170f2-ae87-452f-b1c9-776a6ddfe4b0
plot(0:τ:T, angle.(ψ₁[Nx ÷ 2, :]),
	xlabel = "t",
	title = "arg ψ|<sub>x=0</sub>   (n = $n)",
	xlims = (0, T/2)
)

# ╔═╡ 7d888037-15ac-470e-9f38-dd46dbaaff31
md"""
###### Когерентное состояние
"""

# ╔═╡ 26a8314d-7f00-4940-91d9-0c7157687678
md"""
Когерентное состояние можно выразить через состояния с определённой энергией:
``\displaystyle
|\alpha\rangle = e^{-\alpha^* \alpha / 2} \sum_{n=0}^{\infty} 
\frac{\alpha^n}{\sqrt{n!}} |n\rangle
,\;\;
``
но это неудобно, поскольку придётся иметь дело с бесконечными суммами и вычислять факториалы больших чисел.

Вместо этого будем использовать готовую волновую функцию, которая приведена в [Википедии:](https://en.wikipedia.org/wiki/Coherent_state#The_wavefunction_of_a_coherent_state)
```math
|\alpha\rangle = \pi^{-1/4} \exp 
\biggl(-\frac12 
	\bigl(x - \sqrt2 \operatorname{Re} \alpha \bigr)^2 
 	+ i \sqrt2 \operatorname{Im} \alpha \cdot x
  	- i \operatorname{Im} \alpha \operatorname{Re} \alpha
\biggr)
```

"""

# ╔═╡ 51c8b55f-a9da-4ea0-9aa5-d40501fbb550
# Когерентное состояние (готовая формула)

CoherentState(α, x) = π^(-1/4) * 
	exp(-1/2 * (x - √2 * real(α))^2
		+ im * √2 * imag(α) * x
		- im * real(α) * imag(α)
	);

# ╔═╡ 407ae2c6-c0a9-4f89-afdb-99b2528eefc8
# Когерентное состояние вручную
# (40 членов, при больших α понадобится больше)

CoherentStateBad(α, x) = π^(-1/4) * exp(-abs2(α)/2) *
		Hermite([α^n / factorial(n) / 2^(n/2) for n ∈ big(0) : big(40)])(x) *
		exp(-x^2 / 2);

# ╔═╡ dea1e376-6ad0-424d-be8d-7c7ed429779d
ψ₂ = Array{Complex{Float64}}(undef, Nx, Nt);

# ╔═╡ f07fcc32-fb57-4067-b549-e2c9e86cd946
α = 1 + im

# ╔═╡ b92da645-7afe-418a-bda5-e0cf988f1ff9
begin
	ψ₂[:, 1] .= CoherentState.(α, -L:h:L)
	dsolve!(ψ₂; params...)
end;

# ╔═╡ e8a9b1f2-6168-4a5c-b536-f02d44ab4587
# @bind t₂ Slider(1:Nt, default = 0)
@bind t₂ Clock(0.03, true, false, Nt)

# ╔═╡ d18b415c-916c-43c5-a28e-5a4454b9d8fe
begin
	plot(-L:h:L, 
		[abs.(ψ₂[:, t₂]) angle.(ψ₂[:, t₂])],
		layout = (2,1),
		xlims = (-0.9L, +0.9L),
		ylims = [(-Inf, Inf) (-π, π)],
		labels = ["|ψ|" "arg ψ"],
		ylabel = ["|ψ|" "arg ψ"],
		legend = nothing,
		title = ["|α(t)〉при α(0) = $α и t = $(round((t₂-1)τ, digits = 2))" ""]
	)
end

# ╔═╡ d80e2f0c-0234-41f1-a1b0-386196ed5d53
md"""
Вдали от колокола ``\psi \approx 0``, поэтому там ``\arg \psi`` сильно подвержен ошибкам округления. 

В случае ``\psi = 0`` величина ``\arg \psi`` не определена вовсе.
"""

# ╔═╡ Cell order:
# ╠═3c9376a6-48cb-41a2-bc8a-2706d5b89a52
# ╟─f51b6d70-8dc2-11ec-2f6c-8b52133958a9
# ╟─0ca0b86c-c234-482d-bb39-0be8175a421c
# ╟─672fece0-3f1b-4f1e-a636-2257b16dd2c8
# ╟─db56e6ff-8dd8-456f-8982-7eb11aea7978
# ╟─c2e79179-5bf2-4c79-a293-f3d90d438f95
# ╟─9eac78d8-cb18-4a32-bd64-830f4d041a75
# ╟─82b3aea2-870a-437a-b977-7b94068aaef2
# ╟─9c20ddb6-19a5-4e90-907f-2c55d1a7a5dc
# ╟─139c012b-b1b0-4b93-9739-5ca56230d50e
# ╟─7b45a662-a47d-4b79-8ff9-dcce4254f08e
# ╟─0a56188c-5332-4d97-91fe-a69d901d8925
# ╠═e26bff40-4335-43cb-8332-4ef4203821bb
# ╠═18ea867c-43e5-4d7f-af08-0e73662c393b
# ╠═b7fe7140-8b25-4d7d-86b9-e0436c5af3d9
# ╠═9f70a119-13de-4739-bef9-cf2e1b6dfefe
# ╟─a4a4e044-1956-4efb-8237-b78c67be3174
# ╟─808fc5bb-8eab-4c2e-9772-ef132761ed17
# ╠═666e1ffe-335c-4781-8bc0-a695207f03df
# ╠═97a66676-bd5b-4772-a4ac-3364ea8cb619
# ╠═ba6efd83-bd2a-4870-a350-192d7a213748
# ╠═1ea38beb-bfa5-4f3f-9f58-ce7978f4e087
# ╠═8da9463c-a2b6-4e74-82f8-b8c15e2f5bc8
# ╠═aa2528b7-ea47-44ef-9b18-f9a4f93fd046
# ╟─7566f079-148b-4207-ba77-4f5355b6c394
# ╠═015170f2-ae87-452f-b1c9-776a6ddfe4b0
# ╟─7d888037-15ac-470e-9f38-dd46dbaaff31
# ╟─26a8314d-7f00-4940-91d9-0c7157687678
# ╠═51c8b55f-a9da-4ea0-9aa5-d40501fbb550
# ╠═407ae2c6-c0a9-4f89-afdb-99b2528eefc8
# ╠═dea1e376-6ad0-424d-be8d-7c7ed429779d
# ╠═f07fcc32-fb57-4067-b549-e2c9e86cd946
# ╠═b92da645-7afe-418a-bda5-e0cf988f1ff9
# ╠═e8a9b1f2-6168-4a5c-b536-f02d44ab4587
# ╟─d18b415c-916c-43c5-a28e-5a4454b9d8fe
# ╟─d80e2f0c-0234-41f1-a1b0-386196ed5d53
