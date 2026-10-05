using SciMLBase: ODEProblem, solve
using OrdinaryDiffEq: Tsit5
using Plots

# --------------------------------------------------
# Лабораторная работа №5
# Модель "хищник-жертва"
# Вариант 2
# --------------------------------------------------

# Коэффициенты варианта
a = 0.13
b = 0.33
c = 0.042
d = 0.03

# Начальные численности
# x - хищники
# y - жертвы
x0 = 7.0
y0 = 12.0

# В примере методички используется интервал
# от 0 до 400 с шагом 0.1
t0 = 0.0
tmax = 400.0
dt = 0.1

tspan = (t0, tmax)
save_times = t0:dt:tmax

u0 = [x0, y0]

# Каталог результатов
mkpath("plots/01_predator_prey")

println("="^60)
println("ЛАБОРАТОРНАЯ РАБОТА №5")
println("Модель хищник-жертва")
println("Вариант 2")
println("="^60)

println()
println("Начальные условия:")
println("Хищники x(0) = ", x0)
println("Жертвы   y(0) = ", y0)

println()
println("Коэффициенты:")
println("a = ", a)
println("b = ", b)
println("c = ", c)
println("d = ", d)

# --------------------------------------------------
# Система Лотки-Вольтерры
#
# dx/dt = -0.13x + 0.042xy
# dy/dt =  0.33y - 0.03xy
# --------------------------------------------------

params = (
    a = a,
    b = b,
    c = c,
    d = d
)

function predator_prey!(du, u, p, t)

    x = u[1]
    y = u[2]

    du[1] = -p.a * x + p.c * x * y
    du[2] =  p.b * y - p.d * x * y
end

# --------------------------------------------------
# Стационарное состояние
# --------------------------------------------------
#
# Для ненулевого стационарного состояния:
#
# x* = b / d
# y* = a / c
#

x_stationary = b / d
y_stationary = a / c

println()
println("="^60)
println("СТАЦИОНАРНОЕ СОСТОЯНИЕ")
println("="^60)

println(
    "x* = ",
    round(x_stationary; digits = 4)
)

println(
    "y* = ",
    round(y_stationary; digits = 4)
)

# Проверяем, что производные
# в стационарной точке равны нулю

du_check = zeros(2)

predator_prey!(
    du_check,
    [x_stationary, y_stationary],
    params,
    0.0
)

println()
println("Проверка стационарной точки:")
println(
    "dx/dt = ",
    round(du_check[1]; digits = 10)
)
println(
    "dy/dt = ",
    round(du_check[2]; digits = 10)
)

# --------------------------------------------------
# Численное решение
# --------------------------------------------------

prob = ODEProblem(
    predator_prey!,
    u0,
    tspan,
    params
)

sol = solve(
    prob,
    Tsit5();
    saveat = save_times,
    abstol = 1e-9,
    reltol = 1e-9
)

# Выделяем две популяции
x = [u[1] for u in sol.u]
y = [u[2] for u in sol.u]

println()
println("="^60)
println("ЧИСЛЕННОЕ РЕШЕНИЕ")
println("="^60)

println("Количество временных точек: ", length(sol.t))

println(
    "Диапазон численности хищников: ",
    round(minimum(x); digits = 3),
    " - ",
    round(maximum(x); digits = 3)
)

println(
    "Диапазон численности жертв: ",
    round(minimum(y); digits = 3),
    " - ",
    round(maximum(y); digits = 3)
)

# --------------------------------------------------
# ГРАФИК 1
# Изменение численности хищников
# --------------------------------------------------

p_predators = plot(
    sol.t,
    x,
    xlabel = "Время",
    ylabel = "Численность хищников",
    title = "Изменение численности хищников",
    label = "x(t)",
    lw = 2,
    grid = true
)

savefig(
    p_predators,
    "plots/01_predator_prey/predators_time.png"
)

# --------------------------------------------------
# ГРАФИК 2
# Изменение численности жертв
# --------------------------------------------------

p_prey = plot(
    sol.t,
    y,
    xlabel = "Время",
    ylabel = "Численность жертв",
    title = "Изменение численности жертв",
    label = "y(t)",
    lw = 2,
    grid = true
)

savefig(
    p_prey,
    "plots/01_predator_prey/prey_time.png"
)

# --------------------------------------------------
# ГРАФИК 3
# Зависимость численности хищников
# от численности жертв
# --------------------------------------------------

p_phase = plot(
    y,
    x,
    xlabel = "Численность жертв y",
    ylabel = "Численность хищников x",
    title = "Фазовый портрет системы",
    label = "Траектория",
    lw = 2,
    grid = true
)

# Отмечаем стационарное состояние
scatter!(
    p_phase,
    [y_stationary],
    [x_stationary],
    label = "Стационарное состояние",
    markersize = 6
)

savefig(
    p_phase,
    "plots/01_predator_prey/phase_portrait.png"
)

# --------------------------------------------------
# Все требуемые графики на одном рисунке
# --------------------------------------------------

p_all = plot(
    p_predators,
    p_prey,
    p_phase;
    layout = (1, 3),
    size = (1500, 450)
)

savefig(
    p_all,
    "plots/01_predator_prey/all_results.png"
)

println()
println("="^60)
println("РАСЧЁТ ЗАВЕРШЁН")
println("="^60)

println()
println("Сохранены графики:")
println("plots/01_predator_prey/predators_time.png")
println("plots/01_predator_prey/prey_time.png")
println("plots/01_predator_prey/phase_portrait.png")
println("plots/01_predator_prey/all_results.png")
