# v0.73 文献・資料と採用範囲

2026-11-25確認。原資料へアクセスできた範囲と実装上の近似を区別する。

|資料|確認範囲・採用した知見|実装との対応／限界|
|---|---|---|
|[Teran et al. 2005, Creating and simulating skeletal muscle from the Visible Human Data Set](https://graphics.cs.wisc.edu/Papers/2005/TSBNLF05/)|著者所属機関の論文紹介・要旨。筋束方向、近似非圧縮性、受動／能動成分と接触を含む構成モデル。|受動形状と活動／力を分離する設計。FEM・構成則・30筋接触シミュレーション自体は移植していない。本文を確認済みとは扱わない。|
|[OpenSim Muscle Editor](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53090145/Muscle%20Editor)|公式本文。固定・via・moving・wrap pointと球／楕円体等のwrap object、順序依存性。|付着・経由点・障害物を定義へ分離。現在は基準相対包絡で、OpenSimの接線wrap解法と同等ではない。|
|[Lindenmann et al. 2022, systematic review](https://doi.org/10.1007/s10439-022-03088-8)|公開本文。C1–C2が軸回旋へ主要寄与し、撮像法・対象・姿勢による差がある。|可変設定の関節分配。60%は実測集約値ではなく工学的初期値。個体差や変性状態の推定はしない。|
|[Zhou et al. 2020, C0–T1 in vivo motion](https://pmc.ncbi.nlm.nih.gov/articles/PMC6938025/)|索引化された論文本文・表・考察。PMC直接表示はアクセス制限が出たが、同一論文の検索本文と[PubMed](https://pubmed.ncbi.nlm.nih.gov/31653508/)で照合。上位と下位の寄与は同一ではない。|C1–C2を主要回旋節にする根拠。少人数対象の端点ROMを動的な固定分配と解釈しない。|
|[Wang et al., Computation of Rotation Minimizing Frame](https://www.ag.jku.at/pubs/2007wjzl.pdf)|著者所属機関のPDF本文。最小ねじれフレームと二重反射法。|最小回転による局所横断方向の輸送を採用。論文の高次精度二重反射法全体は未実装。反平行接線は定義した横断軸で処理。|
|[BodyParts3D pinned README](https://github.com/Kevin-Mattheus-Moerman/BodyParts3D/blob/2c83c698788cff57f6317bddcf712e7acfbea7d8/assets/BodyParts3D_data/README_e.html)|元データ説明・リリース3.0表記・使用条件を確認。STL mirrorの部位表、ファイル長、Git SHAも照合。|同一座標系・同一コミット・CC BY-SA 2.1 Japanを維持。近年の別配布版のライセンスへ勝手に置換しない。|
|[Sternocleidomastoid, NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/sites/books/NBK532881/)|本文の読取可能な別URL／検索本文。胸骨柄・鎖骨内側から乳様突起・上項線への付着、多様な頭の構成。|起始を胸骨／鎖骨に分ける構造を用意。個々の頂点の付着は未測定。|
|[Anatomy, Back, Muscles](https://www.ncbi.nlm.nih.gov/books/NBK537074/)|索引本文。頭板状筋と頸板状筋の付着・方向の違い。|別メッシュで保持。頸板状筋の複数頸椎への付着領域同定は残作業。|
|[Wickham & Brown 1998](https://pubmed.ncbi.nlm.nih.gov/9720999/)|原著要旨。三角筋内の複数の解剖学的区画。|前・中・後部を別定義として維持。3区分を実測7区画と同一視せず、羽状角の定量再現は未実装。|

大胸筋の原著については既存[PECTORAL_FAN_V070.md](PECTORAL_FAN_V070.md)の引用を継承。ただし今回のPubMed直接取得は内容制限で本文を再確認できず、寸法や筋束分節を新たな測定値として採用していない。広い起始／収束形態は保持し、個別分節数や遠位腱twistは確定しない。引用の実在と本文内容を全て新規確認できたという報告はしない。
