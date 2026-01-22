#!/bin/bash

# CleanRoom CLI helper script

case $1 in
  "run")
    if [ -z "$2" ]; then
      echo "Usage: ./cr run <param_file.json>"
      exit 1
    fi
    julia --project -e "using CleanroomSolver, CairoMakie; run_simulation(\"$2\")"
    ;;
  "viz")
    STEP=${2:-1000}
    julia --project tools/visualize_step.jl "$STEP"
    ;;
  "test")
    julia --project -e "using Pkg; Pkg.test()"
    ;;
  *)
    echo "CleanRoom CLI helper"
    echo "Usage: ./cr [run|viz|test] [args]"
    echo ""
    echo "Commands:"
    echo "  run <file>    Run simulation with given JSON"
    echo "  viz <step>    Visualize specific step (default: 1000)"
    echo "  test          Run unit tests"
    ;;
esac
