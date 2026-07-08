import React from "react";
import { ExonPrimerPanel } from "./ExonPrimerPanel";

// シーケンス用プライマー設計タブ
// ExonPrimerPanel を Sequencing Mode 固定・文言変更付きで利用する。
export const SequencePrimerPanel: React.FC = () => {
  return (
    <ExonPrimerPanel
      mode="sequence"
      initialSequencingMode
      initialSeqProductMin={600}
      initialSeqProductMax={800}
      initialSeqOverlap={30}
      initialSeqMargin={50}
    />
  );
}

