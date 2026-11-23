# v0.57 胸郭・脊椎の精密化

BodyParts3D / Anatomography (DBCLS), CC BY-SA 2.1 Japan. [Mitsuhashi et al. 2009](https://doi.org/10.1093/nar/gkn613). [STL mirror](https://github.com/Kevin-Mattheus-Moerman/BodyParts3D/tree/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data) / [部位表](https://github.com/Kevin-Mattheus-Moerman/BodyParts3D/blob/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/parts_list_e.txt)。出典リポジトリの部位表とSTLツリーを照合。全48ファイルの存在とバイト数を確認。バイナリはブラウザで取得し、長さと有限座標を検査。

- 右腕と同じ共有座標変換 atlasToWorld を使用。個別メッシュの中心化・回転・拡大を行わない。
- 左右肋骨は実際に別のアトラス部位。左上肢は従来通り鏡映。
- 全ファイル成功時に一括置換。失敗時は全胸郭を旧模式骨格のまま表示し、再試行可。WebGLなしでは軽量表示。
- 胸郭・脊椎は静的。呼吸、関節接触、筋腱付着・滑走、臨床的妥当性は未検証。背部筋はv0.56の近似のまま初期OFF。
- 合成STL・左右操作・表示フィルター・カメラ境界をテスト。実データ実機描画は別途必要。

|ID|File|Group|Bytes|
|---|---|---|---:|
|C1|FMA12519.stl|spine|308684|
|C2|FMA12520.stl|spine|343584|
|C3|FMA12521.stl|spine|263984|
|C4|FMA12522.stl|spine|211284|
|C5|FMA12523.stl|spine|237784|
|C6|FMA12524.stl|spine|317584|
|C7|FMA12525.stl|spine|261184|
|T1|FMA9165.stl|spine|280784|
|T2|FMA9187.stl|spine|418084|
|T3|FMA9209.stl|spine|445584|
|T4|FMA9248.stl|spine|446584|
|T5|FMA9922.stl|spine|491784|
|T6|FMA9945.stl|spine|244884|
|T7|FMA9968.stl|spine|537084|
|T8|FMA9991.stl|spine|260384|
|T9|FMA10014.stl|spine|564384|
|T10|FMA10037.stl|spine|596484|
|T11|FMA10059.stl|spine|282684|
|T12|FMA10081.stl|spine|301984|
|R1|FMA7857.stl|ribs|709784|
|R2|FMA7882.stl|ribs|1184784|
|R3|FMA7909.stl|ribs|1543084|
|R4|FMA7957.stl|ribs|1924484|
|R5|FMA8066.stl|ribs|2081784|
|R6|FMA8175.stl|ribs|2231784|
|R7|FMA8229.stl|ribs|2157384|
|R8|FMA8283.stl|ribs|1832384|
|R9|FMA8364.stl|ribs|1713984|
|R10|FMA8445.stl|ribs|1550584|
|R11|FMA8531.stl|ribs|1121484|
|R12|FMA8533.stl|ribs|663184|
|L1|FMA7987.stl|ribs|773284|
|L2|FMA8012.stl|ribs|1265284|
|L3|FMA8039.stl|ribs|1560684|
|L4|FMA8148.stl|ribs|1896184|
|L5|FMA8093.stl|ribs|2078484|
|L6|FMA8202.stl|ribs|2221784|
|L7|FMA8256.stl|ribs|2203184|
|L8|FMA8310.stl|ribs|1856784|
|L9|FMA8391.stl|ribs|1801084|
|L10|FMA8472.stl|ribs|1587984|
|L11|FMA8532.stl|ribs|1141384|
|L12|FMA8534.stl|ribs|686384|
|manubrium|FMA7486.stl|sternum|594184|
|sternalBody|FMA7487.stl|sternum|995884|
|xiphoid|FMA7488.stl|sternum|204384|
|cartilageR|BP28.stl|cartilage|957884|
|cartilageL|BP24.stl|cartilage|947184|

Total bytes: 48301532
