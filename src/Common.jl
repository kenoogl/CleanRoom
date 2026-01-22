module Common

using FLoops

export DimensionParams, get_backend

"""
    DimensionParams

無次元化パラメータを保持する構造体。
"""
struct DimensionParams
    L0::Float64           # 代表長さ [m]
    U0::Float64           # 代表速度 [m/s]
    nu::Float64           # 動粘性係数 [m²/s]
    Re::Float64           # レイノルズ数 [-]
    T0::Float64           # 代表時間 = L0/U0 [s]
end

"""
    get_backend(par::String)

並列バックエンドを取得する関数。
"thread"の場合は`ThreadedEx()`、それ以外は`SequentialEx()`を返す。
"""
@inline function get_backend(par::String)
    return (par == "thread") ? ThreadedEx() : SequentialEx()
end

end # module Common
