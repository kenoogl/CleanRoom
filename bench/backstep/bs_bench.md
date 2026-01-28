# Benchmark: backward_step



反復法の計時、一度コンパイルのため空打ちしておく

~~~
% time ./cr run verification/backward_step/step.json 
~~~



##### Backward_stepのケースで開発時のベンチ 2026-01-28

| 反復法   | time [s] | コメント    |
| -------- | -------- | ----------- |
| SOR      | 575.18   |             |
| RBSOR    | 467.88   | 2重ループ   |
| RBSOR    | 385.97   | 1重ループ化 |
| SSOR     | 520.59   | 2重ループ   |
| SSOR     | 397.32   | 1重ループ化 |
| RBSSOR   | 279.96   | 2重ループ   |
| CG       | 39.46    | None        |
|          |          | Sor         |
|          |          | Rbsor       |
|          |          | Ssor        |
|          |          | Rbssor      |
| BiCGstab | ?        | None        |
|          | 95.62    | SOR         |
|          | 79.69    | RBSOR       |
|          | 120.76   | SSOR        |
|          | 89.91    | RBSSOR      |

SSORの反復回数は往復を含むので、実際は2倍

反復回数は辞書順のSOR, SSORがベスト



