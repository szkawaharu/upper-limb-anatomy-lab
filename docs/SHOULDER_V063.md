# v0.63.0 肩甲帯・筋変形・起動改善

## 実装
- 右側アトラス座標の前方は −Z。鎖骨は挙上・後退・長軸回旋を独立に近似し、SC–AC距離を保存する。係数は中立から最大挙上までの工学的近似であり、ISB角度や個人実測値ではない。
- 実肋骨点群の高さ別13断面から横断面楕円を推定。肩甲骨内側・下部の点群について中立位からの符号付き半径方向距離変化を最小化。SC/AC接続を保持しつつ4変数を有界探索する。未取得時は基本運動のみ。真の骨接触・貫通回避・靱帯拘束ではない。
- 三角筋・大胸筋の骨頭近傍の半径包絡、前鋸筋の胸郭方向ガード、DQS/線形の局所辺長比較。新モードは中間ウェイト最大30 mmの変位上限、端点の補正は0。旧線形・旧試験モードも選択可。筋の力・体積保存を保証しない。旧線形モードは新肩甲帯姿勢上での線形筋変形比較であり、旧版全体の完全再現ではない。
- 三角筋・大胸筋6部位の閉じた一貫したメッシュのみ符号付き包絡体積を計測。開放面、非多様体、サイズ上限超過は評価不能。自己交差・入れ子の検査は未実施。
- 起動時は簡易メッシュを隠し、全モデル取得まで進捗を表示。未取得の場合は再読み込みか近似モデル表示を選択。固定コミットの胸郭・肩胸背筋STLの一部を、パース成功後に Cache Storage に保持（選定上限48 MiB、1ファイル4 MiB、利用不可時は通常通信）。可視化に必要な処理時間は実測未確認。

## 検証
`node tests/v063-synthetic.mjs`：合成胸郭・肩甲骨点群の63姿勢で有限性・中立位置・鎖骨長・代理誤差非増加を検査。54姿勢で代理誤差低下。経路補正25条件で有限性・端点・変位上限を検査。四面体の閉鎖面・体積と開放面の棄却を検査。JavaScript構文も検査済み。実アトラスの全可動域・WebGLシェーダ実機・iPad Safari・デプロイの確認は別途必要。

## 根拠と限界
- Ludewig et al. (2009), Motion of the shoulder complex during multiplanar humeral elevation. https://pubmed.ncbi.nlm.nih.gov/19181982/ — SC/AC/GHの3次元連動の方向性。個別の係数を転用した実測再現ではない。
- Ludewig et al. (2004), Three-dimensional clavicular motion during arm elevation. https://pubmed.ncbi.nlm.nih.gov/15089027/ — 鎖骨の挙上・後退・後方長軸回旋の根拠。
- Ludewig & Reynolds (2009), The association of scapular kinematics and glenohumeral joint pathologies. https://doi.org/10.2519/JOSPT.2009.2808 — 肩甲帯の運動と臨床的背景。
- Kavan et al. (2008), Geometric Skinning with Approximate Dual Quaternion Blending — GPUスキニングの方法。筋肉の物理モデルや体積保存の証明ではない。
- BodyParts3Dのデータ出典・ライセンスは既存表示を維持。

本版は教育・研究試作。実射の筋張力、障害リスク、正常性の評価には使用しない。
