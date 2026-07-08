import React, { useMemo } from "react";

const clamp = (v: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, v));

const bandY = (size: number, maxSize: number, height: number): number => {
  const s = Math.max(1, size);
  const maxS = Math.max(2, maxSize);
  const minS = 20;
  const logMin = Math.log10(minS);
  const logMax = Math.log10(maxS);
  const t = (Math.log10(s) - logMin) / (logMax - logMin || 1);
  // large fragments stay near top, small fragments go further (down)
  const inv = 1 - clamp(t, 0, 1);
  return clamp(inv * (height - 8) + 2, 2, height - 6);
};

const GelLane: React.FC<{
  fragments: number[];
  maxSize: number;
  height: number;
  width: number;
}> = ({ fragments, maxSize, height, width }) => {
  const bands = useMemo(() => {
    const uniq = Array.from(new Set((fragments ?? []).filter((n) => Number.isFinite(n) && n > 0)));
    uniq.sort((a, b) => b - a);
    return uniq.map((s) => ({
      size: s,
      y: bandY(s, maxSize, height),
    }));
  }, [fragments, height, maxSize]);

  return (
    <div
      style={{
        position: "relative",
        width,
        height,
        borderRadius: 6,
        background: "linear-gradient(180deg, rgba(17,24,39,0.95), rgba(17,24,39,0.75))",
        border: "1px solid rgba(255,255,255,0.08)",
        boxShadow: "inset 0 0 0 1px rgba(0,0,0,0.18)",
        overflow: "hidden",
      }}
      title={bands.map((b) => `${b.size}bp`).join(", ")}
    >
      {bands.map((b) => (
        <div
          key={b.size}
          style={{
            position: "absolute",
            left: 3,
            right: 3,
            top: b.y,
            height: 2,
            borderRadius: 99,
            background: "rgba(255,255,255,0.9)",
            boxShadow: "0 0 10px rgba(255,255,255,0.18)",
          }}
        />
      ))}
    </div>
  );
};

export const GelMini: React.FC<{
  fragmentsRef: number[];
  fragmentsAlt: number[];
  productLenRef?: number;
  productLenAlt?: number;
}> = ({ fragmentsRef, fragmentsAlt, productLenRef, productLenAlt }) => {
  const height = 58;
  const width = 24;
  const maxSize = Math.max(
    100,
    ...(fragmentsRef ?? []),
    ...(fragmentsAlt ?? []),
    productLenRef ?? 0,
    productLenAlt ?? 0,
  );

  if (!fragmentsRef?.length && !fragmentsAlt?.length) return <span>-</span>;

  return (
    <div style={{ display: "flex", gap: 6, alignItems: "center" }}>
      <GelLane fragments={fragmentsRef ?? []} maxSize={maxSize} height={height} width={width} />
      <GelLane fragments={fragmentsAlt ?? []} maxSize={maxSize} height={height} width={width} />
    </div>
  );
};

