using SciMLBase: ODEProblem, solve, ContinuousCallback, terminate!
using OrdinaryDiffEq: Tsit5
using Plots

# --------------------------------------------------
# Лабораторная работа №7
# Эффективность рекламы
# Вариант 2
# --------------------------------------------------

# Общий объём потенциальной аудитории
N = 1000.0

# В начальный момент о товаре знают 2 человека
n0 = 2.0

# Интервал моделирования
t0 = 0.0
tmax = 30.0
dt = 0.001

tspan = (t0, tmax)
save_times = t0:dt:tmax

# Начальное состояние
u0 = [n0]

# Каталог для графиков
mkpath("plots/01_advertising")

println("="^60)
println("ЛАБОРАТОРНАЯ РАБОТА №7")
println("Эффективность рекламы")
println("Вариант 2")
println("="^60)

println()
println("Исходные данные:")
println("N = ", N)
println("n(0) = ", n0)
println("Интервал моделирования: [", t0, "; ", tmax, "]")

# ==================================================
# СЛУЧАЙ 1
#
# dn/dt = (0.65 + 0.0002*n)*(N - n)
# ==================================================

a1 = 0.65
b1 = 0.0002

function advertising1!(du, u, p, t)

    n = u[1]

    du[1] = (
        p.a + p.b * n
    ) * (
        p.N - n
    )

end

params1 = (
    a = a1,
    b = b1,
    N = N
)

prob1 = ODEProblem(
    advertising1!,
    u0,
    tspan,
    params1
)

sol1 = solve(
    prob1,
    Tsit5();
    saveat = save_times,
    abstol = 1e-9,
    reltol = 1e-9
)

n1 = [u[1] for u in sol1.u]

println()
println("="^60)
println("СЛУЧАЙ 1")
println("="^60)

println("alpha1 = ", a1)
println("alpha2 = ", b1)

println(
    "n(30) = ",
    round(n1[end]; digits=4)
)

# График случая 1

p1 = plot(
    sol1.t,
    n1,
    xlabel = "Время",
    ylabel = "Число информированных",
    title = "Случай 1",
    label = "n(t)",
    lw = 3,
    grid = true,
    ylim = (0, N * 1.05)
)

hline!(
    p1,
    [N],
    label = "Размер аудитории N",
    linestyle = :dash
)

savefig(
    p1,
    "plots/01_advertising/case1.png"
)

# ==================================================
# СЛУЧАЙ 2
#
# dn/dt = (0.0003 + 0.9*n)*(N - n)
#
# Дополнительно нужно найти момент,
# когда скорость распространения максимальна.
# ==================================================

a2 = 0.0003
b2 = 0.9

function advertising2!(du, u, p, t)

    n = u[1]

    du[1] = (
        p.a + p.b * n
    ) * (
        p.N - n
    )

end

params2 = (
    a = a2,
    b = b2,
    N = N
)

prob2 = ODEProblem(
    advertising2!,
    u0,
    tspan,
    params2
)

sol2 = solve(
    prob2,
    Tsit5();
    saveat = save_times,
    abstol = 1e-10,
    reltol = 1e-10
)

n2 = [u[1] for u in sol2.u]

# --------------------------------------------------
# Поиск максимальной скорости распространения
# --------------------------------------------------

# Проверяем решение на очень плотной сетке времени

search_times = range(
    t0,
    tmax;
    length = 300001
)

n2_search = [
    sol2(t)[1]
    for t in search_times
]

speed2 = [
    (a2 + b2 * n) * (N - n)
    for n in n2_search
]

max_index = argmax(speed2)

t_max_speed = search_times[max_index]
max_speed = speed2[max_index]
n_at_max_speed = n2_search[max_index]

println()
println("="^60)
println("СЛУЧАЙ 2")
println("="^60)

println("alpha1 = ", a2)
println("alpha2 = ", b2)

println(
    "Момент максимальной скорости: t = ",
    round(t_max_speed; digits=5)
)

println(
    "Максимальная скорость dn/dt = ",
    round(max_speed; digits=4)
)

println(
    "Число информированных в этот момент: n = ",
    round(n_at_max_speed; digits=4)
)

println(
    "n(30) = ",
    round(n2[end]; digits=4)
)

# График распространения для случая 2

p2 = plot(
    sol2.t,
    n2,
    xlabel = "Время",
    ylabel = "Число информированных",
    title = "Случай 2",
    label = "n(t)",
    lw = 3,
    grid = true,
    ylim = (0, N * 1.05)
)

hline!(
    p2,
    [N],
    label = "Размер аудитории N",
    linestyle = :dash
)

savefig(
    p2,
    "plots/01_advertising/case2.png"
)

# График скорости во втором случае

speed2_plot = [
    (a2 + b2 * n) * (N - n)
    for n in n2
]

p2_speed = plot(
    sol2.t,
    speed2_plot,
    xlabel = "Время",
    ylabel = "dn/dt",
    title = "Скорость распространения - случай 2",
    label = "dn/dt",
    lw = 3,
    grid = true,
    xlim = (0, 0.05)
)

vline!(
    p2_speed,
    [t_max_speed],
    label = "Максимальная скорость",
    linestyle = :dash,
    lw = 2
)

savefig(
    p2_speed,
    "plots/01_advertising/case2_speed.png"
)

# ==================================================
# СЛУЧАЙ 3
#
# dn/dt =
# (0.1*sin(2t) + 0.2*cos(3t)*n)*(N - n)
# ==================================================

function advertising3!(du, u, p, t)

    n = u[1]

    alpha1 = 0.1 * sin(2 * t)
    alpha2 = 0.2 * cos(3 * t)

    du[1] = (
        alpha1 + alpha2 * n
    ) * (
        N - n
    )

end

# --------------------------------------------------
# Контроль насыщения аудитории
# --------------------------------------------------
#
# По смыслу модели n не может быть больше N.
# Когда практически вся аудитория уже информирована,
# фиксируем n = N и завершаем интегрирование.
#

saturation_level = N - 1e-6

function saturation_condition(u, t, integrator)
    return u[1] - saturation_level
end

function saturation_affect!(integrator)
    integrator.u[1] = N
    terminate!(integrator)
end

saturation_callback = ContinuousCallback(
    saturation_condition,
    saturation_affect!
)

prob3 = ODEProblem(
    advertising3!,
    u0,
    tspan
)

sol3 = solve(
    prob3,
    Tsit5();
    callback = saturation_callback,
    saveat = save_times,
    abstol = 1e-9,
    reltol = 1e-9
)

# Получаем рассчитанную часть решения
t3 = collect(sol3.t)

n3 = [
    clamp(u[1], 0.0, N)
    for u in sol3.u
]

# Если аудитория достигла насыщения,
# продолжаем график до конца интервала
# постоянным значением n = N.

saturated = abs(sol3.u[end][1] - N) < 1e-5

if saturated

    saturation_time = sol3.t[end]

    if t3[end] < tmax
        push!(t3, tmax)
        push!(n3, N)
    end

end

println()
println("="^60)
println("СЛУЧАЙ 3")
println("="^60)

println("alpha1(t) = 0.1*sin(2t)")
println("alpha2(t) = 0.2*cos(3t)")

if saturated
    println(
        "Аудитория достигла насыщения при t = ",
        round(saturation_time; digits=5)
    )
end

println(
    "Минимальное n(t) = ",
    round(minimum(n3); digits=4)
)

println(
    "Максимальное n(t) = ",
    round(maximum(n3); digits=4)
)

println(
    "n(30) = ",
    round(n3[end]; digits=4)
)

# График случая 3

p3 = plot(
    t3,
    n3,
    xlabel = "Время",
    ylabel = "Число информированных",
    title = "Случай 3",
    label = "n(t)",
    lw = 3,
    grid = true,
    ylim = (0, N * 1.05)
)

hline!(
    p3,
    [N],
    label = "Размер аудитории N",
    linestyle = :dash
)

savefig(
    p3,
    "plots/01_advertising/case3.png"
)

# --------------------------------------------------
# Все три случая
# --------------------------------------------------

p_all = plot(
    p1,
    p2,
    p3;
    layout = (1, 3),
    size = (1600, 480)
)

savefig(
    p_all,
    "plots/01_advertising/all_cases.png"
)

println()
println("="^60)
println("РАСЧЁТ ЗАВЕРШЁН")
println("="^60)

println()
println("Сохранены графики:")
println("plots/01_advertising/case1.png")
println("plots/01_advertising/case2.png")
println("plots/01_advertising/case2_speed.png")
println("plots/01_advertising/case3.png")
println("plots/01_advertising/all_cases.png")
