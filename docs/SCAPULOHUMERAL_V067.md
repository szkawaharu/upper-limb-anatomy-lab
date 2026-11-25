# v0.67.0 — 肩甲骨姿勢レンジ＋GH寄与ガード

## 目的
v0.66で改善した肩甲胸郭の密着性を維持しつつ、高挙上で肩甲骨を回しすぎて胸郭へ合わせる挙動を抑える。さらに、肩甲骨に対する上腕骨の相対挙上が極端にならないよう、GH寄与を広い帯で監査する。

## 設計方針
- v0.66の胸郭面拘束を主ソルバーとして維持。
- v0.67は **v0.66の補正量を後段で監査するガード**。肩甲骨や上腕骨を別の固定軌道へ強制しない。
- 肩甲骨姿勢はアトラス座標の Y-Z-X 分解による proxy 値を使用。ISB臨床角度ではない。
- 既存の基準軌道は、挙上150°付近で上方回旋約50°、後傾proxy約11°となる。これは近年の系統的レビューが示す健常者の広い傾向（150°で上方回旋約48–61°、後傾約11°）と整合するため、そこから広い許容帯を設定。
- 個体差、挙上面差、測定法差が大きいため、単一の正常角度へ固定しない。

## v0.67の処理
1. v0.66が胸郭面に対する候補姿勢を計算。
2. v0.67が未補正姿勢からv0.66候補までを 12.5% 刻みで補間。
3. 各候補で以下を確認。
   - 胸郭幾何コストが未補正姿勢より悪化しない。
   - 肩甲骨 proxy（上方回旋・傾斜・内外旋）が基準軌道から広い許容帯内。
   - GH相対挙上 proxy が挙上角・挙上面に対する広い監査帯内。
   - v0.66による総肩甲骨補正角が挙上角依存の上限内。
4. 条件を満たす最大の補間率を採用。満たさなければ未補正姿勢へ戻す。

## proxy許容帯
基準中心:
- 上方回旋 proxy: humerothoracic elevation / 3
- 傾斜 proxy: 0 → 約13°
- 内外旋 proxy: 挙上面に応じた既存軌道

許容幅:
- 上方回旋: ±(10–15°)
- 傾斜: ±(8–13°)
- 内外旋: ±(14–18°)
- v0.66総補正角: 約10–20°

これらは臨床正常範囲ではなく、既存モデルが文献傾向から大きく逸脱しないための工学的ガード。

## GH監査
- 上腕骨の命令姿勢 qH は変更しない。
- 肩甲骨座標に変換した上腕軸と中立上腕軸のなす角を GH相対挙上 proxy とする。
- 挙上面によって scapulohumeral rhythm が変わることが報告されているため、前方挙上と外転で同じ比率を強制しない。
- 監査帯は意図的に広く設定し、さらに未補正姿勢からのGH proxy変化も制限。
- 生体の上腕骨頭は数mm〜関節窩寸法の一部に相当するtranslationを持つため、本版の固定GH中心は教育用剛体近似のまま。translationそのものは未実装。

## 文献
1. Scapular kinematics variability in individuals with and without rotator cuff-related shoulder pain: systematic review with multilevel meta-regression. PMID 40845626.
   - 健常者の肩甲骨姿勢の大きなばらつきと、挙上に伴う上方回旋・後傾・外旋傾向。150°で上方回旋47.88–61.00°、後傾10.78–11.96°という統合推定を参考にした。
2. Giphart JE, et al. Effect of plane of arm elevation on glenohumeral kinematics: a normative biplane fluoroscopy study. J Bone Joint Surg Am. 2013. PMID 23389787.
   - scapulohumeral rhythmが外転、scaption、前方挙上で異なることを支持。
3. Kuniki M, et al. Shoulder kinematics and muscle synergy during multi-plane humeral elevation and lowering. J Biomech. 2025. PMID 40315759.
   - 挙上面でGH挙上量が異なり、肩甲帯筋活動も関係することを支持。
4. Beyond Euler/Cardan analysis: True glenohumeral axial rotation during arm elevation and rotation. PMID 33989999.
   - Euler/Cardan見かけ角とtrue axial rotationが異なることを示す。v0.67のproxy角を臨床角と呼ばない根拠。
5. Kinematic coupling of the glenohumeral and scapulothoracic joints generates humeral axial rotation. PMID 35367838.
   - STとGHの結合が上腕軸回旋に寄与することを支持。
6. Herbst EC, et al. Determination of inter-individual differences in humeral head centre translations in subjects with healthy shoulders. R Soc Open Sci. 2026. PMID 42841165.
   - 健常者でもHHC translationに個人差・動作差がある。固定GH中心を生体内translationとみなさない根拠。

## 回帰テスト
`node tests/v067-scapulohumeral.mjs`
- 7挙上角 × 5挙上面 × 3軸回旋 = 105条件
- 基準姿勢の誤検出 0
- 有限値
- 胸郭幾何コスト非悪化
- proxy帯内
- GH寄与帯内
- 補正角上限内

開発時の合成検査では105条件中72条件でv0.66補正を部分的に弱め、最大採用補正角は約16°だった。

## 限界
- proxy角はアトラス座標系でありISB角ではない。
- GH中心translation、関節軟骨、接触path、関節包、靱帯、筋力・筋活性は未実装。
- 正常/異常判定や臨床診断には使用しない。
- iPad/iPhone実機スクリーンショットと独立動態画像データでの継続検証が必要。
