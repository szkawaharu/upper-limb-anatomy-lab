# v0.65.0 — 肩甲胸郭面拘束と前鋸筋9スリップ

## 目的
v0.64で残った、挙上時の肩甲骨内側縁・下角の浮きと、前鋸筋が大きな板・トゲ状に変形する問題を対象にする。見た目だけを合わせるのではなく、解剖学的な区域、既存のBodyParts3D表面、胸郭・肩甲骨の相対関係を使う。

## 肩甲胸郭
- 肩甲骨肋骨面のサンプルを使い、実肋骨点群から作った高さ別胸郭断面に対する中立位相対距離を評価する。
- SC–AC接続と鎖骨長を維持する既存5変数探索を継続。
- v0.65では平均二乗誤差だけでなく、局所的な最大離開と内方移動も罰則化する。8 mm（外方）・4 mm（内方）からの追加罰則は数値安定化のための工学的閾値であり、正常範囲や臨床的クリアランスではない。
- 真の胸郭面接触、肩甲下筋・前鋸筋厚、滑液包、軟部組織圧縮、鎖骨/AC/SCの関節力は未モデル化。

## 前鋸筋
- 文献では通常、第1〜8または第1〜9肋骨から複数のserration/slipとして起始し、肩甲骨内側縁の肋骨面に停止する。上部・中部・下部の3区域が記載される。
- v0.65では第1〜9肋骨を9本の連続スリップ座標として扱う。各スリップの肋骨側起始は、BodyParts3D前鋸筋の起始側表面と各実肋骨サンプルとの最近接候補から決定する。
- 肩甲骨側は、既存の内側縁ランドマーク経路を上部（1–2）、中部（3–6）、下部（6–9）の区域に沿って連続補間する。個体の実測停止面ではない。
- 各スリップは4点の曲線で表し、中間2点は中立位の実筋メッシュから推定。隣接スリップ間を連続補間して継ぎ目を避ける。
- 変形後は既存v0.64方式を基準に局所辺長と胸郭包絡代理を比較し、悪化する候補は採用しない。v0.65 GPUシェーダが使えない端末ではv0.64へ戻す。
- 9本は「解剖学的区域の表示単位」であり、実測された筋束数・運動単位・独立した力ベクトルを意味しない。深浅subslipも未分離。

## 拡張性
serratus-v065.js を独立モジュール化。将来、肋骨別の停止面、深浅subslip、長胸神経支配、実測筋束方向を追加しても、共通筋変形エンジンを壊さず専用定義を更新できる。

## 検証
- `node tests/v065-serratus.mjs`: 合成表面で9スリップ、連続座標、中立形状、決定性、候補非悪化を検査。
- `node tests/v065-shoulder.mjs`: 合成胸郭で105姿勢、有限性、SC–AC長、補正目的関数の非悪化を検査。
- 実際のiOS SafariにおけるGLSLのv0.65経路、全17筋の実形状、前鋸筋の肋骨別停止、肩甲骨の実接触距離は実機画像で継続検証する。

## 主な解剖学・生体力学資料
1. Seth A, et al. A Biomechanical Model of the Scapulothoracic Joint to Accurately Capture Scapular Kinematics during Shoulder Movements. PLoS ONE. 2016;11:e0141028. https://doi.org/10.1371/journal.pone.0141028
   - 胸郭表面に対する肩甲骨の並進・回転を拘束する考え方。本実装は同論文モデルの再実装ではない。
2. Ludewig PM, et al. Motion of the Shoulder Complex During Multiplanar Humeral Elevation. J Bone Joint Surg Am. 2009. PMID 19181982.
   - SC、AC、ST、GHの連動を扱う基礎。個人別角度を直接転用していない。
3. Lung K, St Lucia K, Lui F. Anatomy, Thorax, Serratus Anterior Muscles. StatPearls. PMID 30285352.
   - 第1〜8/9肋骨から肩甲骨内側縁への扇状起始・停止。
4. Serratus anterior palsy secondary to long thoracic nerve dysfunction. 2025 review.
   - 上部1–2肋骨、中部3–6肋骨、下部6以下という区域と肩甲骨内側縁・下角への対応。
5. Manktelow RT, et al. Detailed neurovascular anatomy of the serratus anterior muscle: implications for a functional muscle flap with multiple independent force vectors. Plast Reconstr Surg. 2004. PMID 15220562.
   - 下位スリップの独立性を支持。ただし本モデルは血管・神経単位を再現しない。
6. Hamada J, et al. Anatomical characteristics of the upper serratus anterior: cadaver dissection. Clin Anat. 2003. PMID 12968858.
   - 上部前鋸筋の個体差と独立した形態。9本を固定した普遍的解剖と見なさない根拠。

## 医学教育上の位置付け
v0.65は研究試作。医学部教育にそのまま採用できる検証水準ではない。特に、前鋸筋の停止面・subslip、肩甲骨と胸郭の軟部組織厚、個体差、動態画像との照合、専門家レビューが必要。
