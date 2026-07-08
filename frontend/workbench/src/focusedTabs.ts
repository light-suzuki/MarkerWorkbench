import React from "react";

const ExonPrimerPanel = React.lazy(async () => ({ default: (await import("./components/ExonPrimerPanel")).ExonPrimerPanel }));
const SequencePrimerPanel = React.lazy(async () => ({ default: (await import("./components/SequencePrimerPanel")).SequencePrimerPanel }));
const CapsPrimerPanel = React.lazy(async () => ({ default: (await import("./components/CapsPrimerPanel")).CapsPrimerPanel }));

export const productName = "Marker Workbench";
export const focusedTabs = [
  { id: "exon", labelJa: "CDS/エキソン", labelEn: "CDS/exon primers", descriptionJa: "CDS・エキソン増幅候補を設計", descriptionEn: "Design CDS and exon amplification candidates", color: "#ea580c", Component: ExonPrimerPanel },
  { id: "seq_primers", labelJa: "シーケンスプライマー", labelEn: "Sequencing primers", descriptionJa: "シーケンス用増幅候補を設計", descriptionEn: "Design sequencing amplicons", color: "#0d9488", Component: SequencePrimerPanel },
  { id: "caps", labelJa: "CAPS", labelEn: "CAPS markers", descriptionJa: "CAPSマーカー候補を生成", descriptionEn: "Generate CAPS marker candidates", color: "#b45309", Component: CapsPrimerPanel },
] as const;

