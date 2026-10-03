#!/usr/bin/env julia

using Pkg

Pkg.activate(".")

packages = [
    "DrWatson",
    "DifferentialEquations",
    "Plots",
    "DataFrames",
    "CSV",
    "JLD2",
    "Literate",
    "IJulia",
    "BenchmarkTools",
    "Quarto"
]

println("Установка базовых пакетов...")

Pkg.add(packages)

println()
println("Все пакеты установлены.")
println("Для проверки: using DrWatson, DifferentialEquations, Plots")
