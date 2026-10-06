using SciMLBase: ODEProblem, solve
using OrdinaryDiffEq: Tsit5
using Plots

# --------------------------------------------------
# Лабораторная работа №8
# Модель конкуренции двух фирм
# Вариант 2
# --------------------------------------------------

# Начальные значения оборотных средств, млн ед.
M10 = 1.9
M20 = 2.1

# Параметры варианта
# p_cr, p1, p2 и N заданы в тысячах единиц
p_cr = 25.0
N = 50.0
q = 1.0

tau1 = 30.0
tau2 = 20.0

# Себестоимость продукции фирм
p1 = 7.0
p2 = 10.0

# Дополнительный коэффициент
# социально-психологического воздействия
social = 0.0012

# --------------------------------------------------
# Коэффициенты модели
# --------------------------------------------------

a1 = p_cr / (
    tau1^2 * p1^2 * N * q
)

a2 = p_cr / (
    tau2^2 * p2^2 * N * q
)

b = p_cr / (
    tau1^2 * p1^2 *
    tau2^2 * p2^2 *
    N * q
)

c1 = (
    p_cr - p1
) / (
    tau1 * p1
)

c2 = (
    p_cr - p2
) / (
    tau2 * p2
)

# --------------------------------------------------
# Временной интервал
# theta = t / c1
# --------------------------------------------------

theta0 = 0.0
theta_max = 30.0
dtheta = 0.01

theta_span = (
    theta0,
    theta_max
)

save_times = theta0:dtheta:theta_max

# Начальный вектор
u0 = [
    M10,
    M20
]

# Каталог графиков
mkpath("plots/01_firm_competition")

println("="^65)
println("ЛАБОРАТОРНАЯ РАБОТА №8")
println("Модель конкуренции двух фирм")
println("Вариант 2")
println("="^65)

println()
println("Исходные данные:")
println("M1(0) = ", M10, " млн ед.")
println("M2(0) = ", M20, " млн ед.")
println("p_cr = ", p_cr)
println("N = ", N)
println("q = ", q)
println("tau1 = ", tau1)
println("tau2 = ", tau2)
println("p1 = ", p1)
println("p2 = ", p2)

println()
println("Рассчитанные коэффициенты:")
println("a1 = ", a1)
println("a2 = ", a2)
println("b  = ", b)
println("c1 = ", c1)
println("c2 = ", c2)

# ==================================================
# СЛУЧАЙ 1
#
# Конкуренция только рыночными методами
# ==================================================

function competition_case1!(du, u, p, theta)

    M1 = u[1]
    M2 = u[2]

    du[1] =
        M1 -
        (p.b / p.c1) * M1 * M2 -
        (p.a1 / p.c1) * M1^2

    du[2] =
        (p.c2 / p.c1) * M2 -
        (p.b / p.c1) * M1 * M2 -
        (p.a2 / p.c1) * M2^2
end

params1 = (
    a1 = a1,
    a2 = a2,
    b = b,
    c1 = c1,
    c2 = c2
)

prob1 = ODEProblem(
    competition_case1!,
    u0,
    theta_span,
    params1
)

sol1 = solve(
    prob1,
    Tsit5();
    saveat = save_times,
    abstol = 1e-9,
    reltol = 1e-9
)

M1_case1 = [
    u[1] for u in sol1.u
]

M2_case1 = [
    u[2] for u in sol1.u
]

println()
println("="^65)
println("СЛУЧАЙ 1")
println("Только рыночная конкуренция")
println("="^65)

println(
    "M1(30) = ",
    round(M1_case1[end]; digits=4),
    " млн ед."
)

println(
    "M2(30) = ",
    round(M2_case1[end]; digits=4),
    " млн ед."
)

# График первого случая

p_case1 = plot(
    sol1.t,
    M1_case1,
    xlabel = "Безразмерное время θ",
    ylabel = "Оборотные средства, млн ед.",
    title = "Случай 1: рыночная конкуренция",
    label = "Фирма 1",
    lw = 3,
    grid = true
)

plot!(
    p_case1,
    sol1.t,
    M2_case1,
    label = "Фирма 2",
    lw = 3
)

savefig(
    p_case1,
    "plots/01_firm_competition/case1.png"
)

# ==================================================
# СЛУЧАЙ 2
#
# Помимо рыночной конкуренции учитывается
# социально-психологический фактор.
#
# В первом уравнении коэффициент при M1*M2
# увеличивается на 0.0012.
# ==================================================

function competition_case2!(du, u, p, theta)

    M1 = u[1]
    M2 = u[2]

    du[1] =
        M1 -
        (
            p.b / p.c1 +
            p.social
        ) * M1 * M2 -
        (p.a1 / p.c1) * M1^2

    du[2] =
        (p.c2 / p.c1) * M2 -
        (p.b / p.c1) * M1 * M2 -
        (p.a2 / p.c1) * M2^2
end

params2 = (
    a1 = a1,
    a2 = a2,
    b = b,
    c1 = c1,
    c2 = c2,
    social = social
)

prob2 = ODEProblem(
    competition_case2!,
    u0,
    theta_span,
    params2
)

sol2 = solve(
    prob2,
    Tsit5();
    saveat = save_times,
    abstol = 1e-9,
    reltol = 1e-9
)

M1_case2 = [
    u[1] for u in sol2.u
]

M2_case2 = [
    u[2] for u in sol2.u
]

# Максимальное значение фирмы 1 во втором случае

peak_index = argmax(M1_case2)

M1_peak = M1_case2[peak_index]
theta_peak = sol2.t[peak_index]

println()
println("="^65)
println("СЛУЧАЙ 2")
println("С учётом социально-психологического фактора")
println("="^65)

println("Дополнительный коэффициент = ", social)

println(
    "Максимум M1 = ",
    round(M1_peak; digits=4),
    " млн ед."
)

println(
    "Момент максимума фирмы 1: theta = ",
    round(theta_peak; digits=4)
)

println(
    "M1(30) = ",
    round(M1_case2[end]; digits=6),
    " млн ед."
)

println(
    "M2(30) = ",
    round(M2_case2[end]; digits=4),
    " млн ед."
)

# График второго случая

p_case2 = plot(
    sol2.t,
    M1_case2,
    xlabel = "Безразмерное время θ",
    ylabel = "Оборотные средства, млн ед.",
    title = "Случай 2: социальный фактор",
    label = "Фирма 1",
    lw = 3,
    grid = true
)

plot!(
    p_case2,
    sol2.t,
    M2_case2,
    label = "Фирма 2",
    lw = 3
)

savefig(
    p_case2,
    "plots/01_firm_competition/case2.png"
)

# --------------------------------------------------
# Сравнение двух случаев
# --------------------------------------------------

p_comparison = plot(
    p_case1,
    p_case2;
    layout = (1, 2),
    size = (1400, 520)
)

savefig(
    p_comparison,
    "plots/01_firm_competition/comparison.png"
)

println()
println("="^65)
println("РАСЧЁТ ЗАВЕРШЁН")
println("="^65)

println()
println("Сохранены графики:")
println("plots/01_firm_competition/case1.png")
println("plots/01_firm_competition/case2.png")
println("plots/01_firm_competition/comparison.png")
