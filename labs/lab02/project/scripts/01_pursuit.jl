using SciMLBase: ODEProblem, solve
using OrdinaryDiffEq: Tsit5
using Plots
using DataFrames
using CSV

# --------------------------------------------------
# Лабораторная работа №2
# Вариант 2
# --------------------------------------------------

# Начальное расстояние между лодкой и катером
k = 12.0

# Катер движется в 4 раза быстрее лодки
n = 4.0

# Направление движения лодки
# В примере методички используется 3π/4
phi = 3pi / 4

println("ЛАБОРАТОРНАЯ РАБОТА №2")
println("Вариант 2")
println()
println("k = ", k, " км")
println("Скорость катера в ", n, " раза больше скорости лодки")

# --------------------------------------------------
# ЗАДАНИЕ 1
# Начальные условия для двух случаев
# --------------------------------------------------

# В первом случае:
# x/v = (k-x)/(nv)
# поэтому x = k/(n+1)

r01 = k / (n + 1)

# Во втором случае:
# x/v = (k+x)/(nv)
# поэтому x = k/(n-1)

r02 = k / (n - 1)

println()
println("ЗАДАНИЕ 1")
println("Первый случай:")
println("theta0 = 0")
println("r0 = ", r01, " км")

println()
println("Второй случай:")
println("theta0 = -pi")
println("r0 = ", r02, " км")

# Скорость катера раскладывается
# на радиальную и тангенциальную части.
#
# dr/dtheta = r / sqrt(n^2 - 1)

a = sqrt(n^2 - 1)

println()
println("Уравнение движения катера:")
println("dr/dtheta = r / sqrt(15)")

# --------------------------------------------------
# Дифференциальное уравнение
# --------------------------------------------------

function pursuit!(du, u, p, theta)
    du[1] = u[1] / p
end

# --------------------------------------------------
# СЛУЧАЙ 1
# theta0 = 0
# r0 = 2.4
# --------------------------------------------------

prob1 = ODEProblem(
    pursuit!,
    [r01],
    (0.0, phi),
    a
)

sol1 = solve(
    prob1,
    Tsit5();
    saveat=0.01
)

theta1 = sol1.t
r1 = first.(sol1.u)

# Перевод из полярных координат
# в декартовы для построения графика

x1 = r1 .* cos.(theta1)
y1 = r1 .* sin.(theta1)

# --------------------------------------------------
# СЛУЧАЙ 2
# theta0 = -pi
# r0 = 4
# --------------------------------------------------

prob2 = ODEProblem(
    pursuit!,
    [r02],
    (-pi, phi),
    a
)

sol2 = solve(
    prob2,
    Tsit5();
    saveat=0.01
)

theta2 = sol2.t
r2 = first.(sol2.u)

x2 = r2 .* cos.(theta2)
y2 = r2 .* sin.(theta2)

# --------------------------------------------------
# ЗАДАНИЕ 3
# Точки пересечения катера и лодки
# --------------------------------------------------

r_meet1 = sol1(phi)[1]
r_meet2 = sol2(phi)[1]

x_meet1 = r_meet1 * cos(phi)
y_meet1 = r_meet1 * sin(phi)

x_meet2 = r_meet2 * cos(phi)
y_meet2 = r_meet2 * sin(phi)

println()
println("ЗАДАНИЕ 3")
println("Точки пересечения:")

println()
println("Случай 1:")
println(
    "x = ", round(x_meet1; digits=4),
    ", y = ", round(y_meet1; digits=4)
)

println()
println("Случай 2:")
println(
    "x = ", round(x_meet2; digits=4),
    ", y = ", round(y_meet2; digits=4)
)

# --------------------------------------------------
# Сохраняем численные результаты
# --------------------------------------------------

results = DataFrame(
    case = ["Случай 1", "Случай 2"],
    r0 = [r01, r02],
    x_meet = [x_meet1, x_meet2],
    y_meet = [y_meet1, y_meet2]
)

mkpath("data/01_pursuit")
CSV.write(
    "data/01_pursuit/intersections.csv",
    results
)

# --------------------------------------------------
# ЗАДАНИЕ 2
# График первого случая
# --------------------------------------------------

p1 = plot(
    xlabel = "x, км",
    ylabel = "y, км",
    title = "Задача о погоне - случай 1",
    aspect_ratio = :equal,
    legend = :topleft,
    grid = true
)

# Катер сначала движется по прямой
plot!(
    p1,
    [k, r01],
    [0.0, 0.0],
    label = "Катер: прямой участок",
    lw = 2
)

# Затем катер движется по криволинейной траектории
plot!(
    p1,
    x1,
    y1,
    label = "Катер: траектория погони",
    lw = 3
)

# Лодка движется прямолинейно
plot!(
    p1,
    [0.0, x_meet1],
    [0.0, y_meet1],
    label = "Лодка",
    lw = 2,
    linestyle = :dash
)

# Точка встречи
scatter!(
    p1,
    [x_meet1],
    [y_meet1],
    label = "Точка встречи",
    markersize = 6
)

# --------------------------------------------------
# График второго случая
# --------------------------------------------------

p2 = plot(
    xlabel = "x, км",
    ylabel = "y, км",
    title = "Задача о погоне - случай 2",
    aspect_ratio = :equal,
    legend = :topleft,
    grid = true
)

# Во втором случае катер проходит
# на противоположную сторону от полюса
plot!(
    p2,
    [k, -r02],
    [0.0, 0.0],
    label = "Катер: прямой участок",
    lw = 2
)

plot!(
    p2,
    x2,
    y2,
    label = "Катер: траектория погони",
    lw = 3
)

plot!(
    p2,
    [0.0, x_meet2],
    [0.0, y_meet2],
    label = "Лодка",
    lw = 2,
    linestyle = :dash
)

scatter!(
    p2,
    [x_meet2],
    [y_meet2],
    label = "Точка встречи",
    markersize = 6
)

# --------------------------------------------------
# Сохраняем графики
# --------------------------------------------------

mkpath("plots/01_pursuit")

savefig(
    p1,
    "plots/01_pursuit/case1.png"
)

savefig(
    p2,
    "plots/01_pursuit/case2.png"
)

p_all = plot(
    p1,
    p2;
    layout = (1, 2),
    size = (1200, 500)
)

savefig(
    p_all,
    "plots/01_pursuit/both_cases.png"
)

println()
println("Графики сохранены:")
println("plots/01_pursuit/case1.png")
println("plots/01_pursuit/case2.png")
println("plots/01_pursuit/both_cases.png")
