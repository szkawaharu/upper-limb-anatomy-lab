# v0.57.1 — 胸郭接続の補完

出典：BodyParts3D / Anatomography (DBCLS), CC BY-SA 2.1 Japan. [Mitsuhashi et al. 2009](https://doi.org/10.1093/nar/gkn613). [固定した上流部位一覧](https://github.com/Kevin-Mattheus-Moerman/BodyParts3D/blob/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/parts_list_e.txt) および同一コミットのSTLツリーを照合。

v0.57ではBP24/28のみ読み込み、第1〜7肋軟骨のFMA個別ファイルが未読込だった。旧BP部位の正確な構成範囲はメタデータのみでは確定できないため残し、個別表示・実メッシュ三角形交差監査を追加した。交差候補は自動削除しない（正常な接触を削る危険）。剣状突起は従来のFMA7488を保持し、骨化状態を推測して追加軟骨を生成しない。

椎間板はC2〜C7、T1〜T11に対応する17個の下位椎間板。C1–C2間には椎間板がなく、T12–L1は腰椎未実装につき除外。データの元座標を保持。

基本48部位・追加肋軟骨14・椎間板17をグループ単位で置換。追加グループ取得に失敗しても基本胸郭は保持し、不足グループのみ再取得。

検査：旧左右軟骨と新しい左右第1〜7肋軟骨の三角形交差候補、7番軟骨と第7肋骨・旧軟骨・胸骨体・剣状突起のAABB重なり、左右SC/AC推定中心と中立時メッシュ頂点距離。包含・共面重複・関節面接触・人体精度は未検証。実データの結果は端末で検査ボタンを押して確認する。関節中心を未検証のまま移動しない。

|追加ID|ファイル|バイト|
|---|---|---:|
|cartR1|FMA7875.stl|117184|
|cartR2|FMA7886.stl|127084|
|cartR3|FMA7913.stl|134784|
|cartR4|FMA7976.stl|174384|
|cartR5|FMA8070.stl|241584|
|cartR6|FMA8194.stl|303984|
|cartR7|FMA8248.stl|524684|
|cartL1|FMA8005.stl|110284|
|cartL2|FMA8031.stl|117584|
|cartL3|FMA8058.stl|121484|
|cartL4|FMA8167.stl|161084|
|cartL5|FMA8112.stl|226884|
|cartL6|FMA8221.stl|298384|
|cartL7|FMA8275.stl|495184|
|discC2|FMA25058.stl|150784|
|discC3|FMA13896.stl|151684|
|discC4|FMA13897.stl|175884|
|discC5|FMA13898.stl|172284|
|discC6|FMA13899.stl|161984|
|discC7|FMA13900.stl|195684|
|discT1|FMA10458.stl|214084|
|discT2|FMA13495.stl|248584|
|discT3|FMA13500.stl|258784|
|discT4|FMA13501.stl|315484|
|discT5|FMA13502.stl|366884|
|discT6|FMA13503.stl|452984|
|discT7|FMA13504.stl|519384|
|discT8|FMA13505.stl|510184|
|discT9|FMA13506.stl|252084|
|discT10|FMA13507.stl|576784|
|discT11|FMA13508.stl|714084|

全79ファイル合計 56893736 bytes。合成STLと取得失敗を含む統合・既存回帰試験。実機描画未確認。
