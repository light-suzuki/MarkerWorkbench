import React, { useEffect, useMemo, useState } from "react";
import { bioapiClient } from "../api/bioapiClient";
import type { CapsDesignResponse, CapsMarkerRow } from "../types/caps";
import { downloadMarkdown, openPrintViewForMarkdown } from "../utils/exportReport";
import { downloadXlsx } from "../utils/exportXlsx";
import type { XlsxSheet } from "../utils/exportXlsx";
import { ensemblGeneUrl, isLocalOnlyDb, navigatorGeneUrl } from "../utils/ensembl";
import {
  CUSTOM_DB_VALUE,
  FALLBACK_LOCAL_DB_OPTIONS,
  query_TO_ref_VIRTUAL_DB_VALUE,
  labelForDbPath,
  loadPreferredLocalDbPaths,
  normalizeLocalDbValue,
  useLocalBlastDbOptions,
  withCustomDbOption,
} from "../utils/localBlastDbs";
import type { BlastDbChromosome } from "../types/blast";
import type { JobInfo } from "../types/jobs";
import { GelMini } from "./GelMini";
import { useToast } from "./ToastProvider";
import { pollJobUntilDone } from "../utils/jobPolling";
import { JobProgressCard } from "./JobProgressCard";
import { EnsemblLinksInline } from "./EnsemblLinksInline";

const DEFAULT_ENZYMES_TEXT = [
  "EcoRI",
  "HindIII",
  "BamHI",
  "PstI",
  "XbaI",
  "SpeI",
  "XhoI",
  "SalI",
  "SacI",
  "KpnI",
  "ApaI",
  "NotI",
  "EcoRV",
  "BglII",
  "NcoI",
  "NdeI",
  "HaeIII",
  "MspI",
  "RsaI",
  "AluI",
  "HinfI",
  "DdeI",
  "TaqI",
  "Sau3AI",
  "MboI",
  "AvaII",
  "BstUI",
  "BfaI",
  "NlaIII",
  "MseI",
  "PvuII",
  "SmaI",
  "SspI",
  "DraI",
  "ClaI",
  "EcoO109I",
  "AflII",
  "AgeI",
  "BsrGI",
  "BsiWI",
  "BstEII",
  "NheI",
  "SbfI",
].join("\n");

const labelForDb = (path: string) =>
  labelForDbPath(path, null);

const normalizeEnzymeList = (text: string): string[] =>
  text
    .split(/[\s,]+/g)
    .map((t) => t.trim())
    .filter(Boolean);

const formatFragments = (frags: number[]) =>
  frags && frags.length ? frags.join("+") : "-";

const formatRange = (start: number, end: number) =>
  `${start.toLocaleString()}–${end.toLocaleString()}`;

const normalizeChr = (raw: string): string => {
  const t = (raw ?? "").trim();
  if (!t) return "";
  const m = /^(?:chr)?0*([0-9]+)$/i.exec(t);
  if (m) return `chr${Number(m[1])}`;
  if (/^chr/i.test(t)) return `chr${t.slice(3)}`;
  return t;
};

const resolveEntryFromChrom = (raw: string, chroms: BlastDbChromosome[] | null): string => {
  const t = (raw ?? "").trim();
  if (!t) return "";
  const n = normalizeChr(t);
  const hit = (chroms ?? []).find((c) => normalizeChr(c.chrom) === n);
  return hit?.entry || t;
};

const buildTsv = (res: CapsDesignResponse): string => {
  const dbs = Array.from(
    new Set(res.markers.flatMap((m) => (m.blast ?? []).map((b) => b.db))),
  );
  const header = [
    "#",
    "enzyme",
    "gene",
    "ref_range",
    "alt_range",
    "alt_strand",
    "product_ref",
    "product_alt",
    "fragments_ref",
    "fragments_alt",
    "mismatch",
    "primer_left",
    "primer_right",
    ...dbs.flatMap((db) => [
      `${db}_amplicons`,
      `${db}_quality`,
      `${db}_top`,
      `${db}_gene`,
    ]),
  ].join("\t");

  const lines = res.markers.map((m) => {
    const row = [
      m.index,
      m.enzyme,
      m.gene_label ?? "",
      formatRange(m.ref_product_start, m.ref_product_end),
      formatRange(m.alt_product_start, m.alt_product_end),
      m.alt_strand,
      m.product_len_ref,
      m.product_len_alt,
      formatFragments(m.fragments_ref),
      formatFragments(m.fragments_alt),
      m.mismatch_count,
      m.primer_left,
      m.primer_right,
    ].map(String);

    const map = new Map((m.blast ?? []).map((b) => [b.db, b]));
    dbs.forEach((db) => {
      const b = map.get(db);
      const top =
        b?.top_subject && b.top_start != null && b.top_end != null
          ? `${b.top_subject}:${formatRange(b.top_start, b.top_end)}`
          : "";
      row.push(String(b?.amplicon_count ?? ""));
      row.push(String(b?.quality ?? ""));
      row.push(top);
      row.push(String(b?.gene_label ?? ""));
    });
    return row.join("\t");
  });

  return [header, ...lines].join("\n");
};

const buildMarkdown = (res: CapsDesignResponse): string => {
  const lines: string[] = [];
  const dt = new Date();
  lines.push("# CAPS プライマー作成レポート");
  lines.push("");
  lines.push(`- 作成時刻: ${dt.toLocaleString()}`);
  lines.push(`- 参照: \`${labelForDb(res.ref_db)}\` / \`${res.ref_entry}\` / ${res.ref_start}–${res.ref_end} (${res.ref_length} bp)`);
  lines.push(
    `- 比較: \`${labelForDb(res.alt_db)}\` / \`${res.alt_entry}\` / ${res.alt_start}–${res.alt_end} (${res.alt_length} bp) / strand=${res.alt_strand} / mapped_by_blast=${res.mapped_by_blast}`,
  );
  lines.push(`- Primer3 生成ペア数: ${res.primer_pairs_generated}`);
  lines.push(`- 返却マーカー数: ${res.markers.length}`);
  if (res.warnings?.length) {
    lines.push("");
    lines.push("## Warnings");
    lines.push("");
    res.warnings.forEach((w) => lines.push(`- ${w}`));
  }
  lines.push("");
  lines.push("## Markers");
  lines.push("");
  lines.push("| # | Enzyme | Gene | Ref range | Alt range | Ref frag | Alt frag | mismatch |");
  lines.push("| ---: | --- | --- | --- | --- | --- | --- | ---: |");
  res.markers.slice(0, 200).forEach((m) => {
    const url = ensemblGeneUrl(m.gene_label ?? null);
    // TODO: markdown report logic update for Navigator (simplified for now to stick to Ensembl for plain text report or update helper?)
    // For Markdown, we might want to keep simple links. Let's rely on the assumption that if it's local only, we might just put the ID or a Navigator link if possible.
    // For now, let's leave md simple or use the new helper if we want.
    // Actually, let's use the helper to get the RIGHT url.
    // But we need to know if it's local db. Caps result has blast details but here we just have gene_label.
    // We can infer from the DB context or just check the label pattern?
    // gene_label alone might be GENE05G... which is Ensembl.
    // If it's something else, maybe Navigator.
    // For safety in Markdown, let's stick to ensemblGeneUrl for GENE IDs and maybe nothing for others unless we want to embed localhost links which might not work for sharing.
    // Let's leave Markdown as is (Ensembl only) to valid external links, or maybe update if user requests.
    // User request was "ref/ALT except -> Navigator".
    // So if I have a gene label that is NOT Ensembl, I should link to Navigator?
    // Markdown is static text, localhost link is fine for personal use.
    // Let's update it.
    /*
    const localOnly = isLocalOnlyDb(...); // We don't have DB context easily here per gene...
    Actually `m` comes from `res.markers`. `res` has `ref_db`.
    But genes come from BLAST hits against MANY DBs.
    Wait, `m.gene_label` is the primary gene label.
    */

    const gene = m.gene_label ? (url ? `[${m.gene_label}](${url})` : m.gene_label) : "-";
    lines.push(
      `| ${m.index} | ${m.enzyme} | ${gene} | ${formatRange(m.ref_product_start, m.ref_product_end)} | ${formatRange(m.alt_product_start, m.alt_product_end)} (${m.alt_strand}) | ${formatFragments(m.fragments_ref)} | ${formatFragments(m.fragments_alt)} | ${m.mismatch_count} |`,
    );
  });
  if (res.markers.length > 200) {
    lines.push("");
    lines.push(`※ 表は先頭 200 件まで（全 ${res.markers.length} 件）`);
  }
  lines.push("");
  return lines.join("\n");
};

const buildXlsxSheets = (res: CapsDesignResponse): XlsxSheet[] => {
  const dbs = Array.from(new Set(res.markers.flatMap((m) => (m.blast ?? []).map((b) => b.db))));
  const header: Array<string> = [
    "#",
    "enzyme",
    "gene",
    "ref_start",
    "ref_end",
    "alt_start",
    "alt_end",
    "alt_strand",
    "product_len_ref",
    "product_len_alt",
    "fragments_ref",
    "fragments_alt",
    "mismatch",
    "primer_left",
    "primer_right",
    ...dbs.flatMap((db) => [
      `${db}_amplicons`,
      `${db}_quality`,
      `${db}_top_subject`,
      `${db}_top_start`,
      `${db}_top_end`,
      `${db}_gene`,
    ]),
  ];

  const rows = res.markers.map((m) => {
    const row: Array<string | number | null> = [
      m.index,
      m.enzyme,
      m.gene_label ?? "",
      m.ref_product_start,
      m.ref_product_end,
      m.alt_product_start,
      m.alt_product_end,
      m.alt_strand,
      m.product_len_ref,
      m.product_len_alt,
      formatFragments(m.fragments_ref),
      formatFragments(m.fragments_alt),
      m.mismatch_count,
      m.primer_left,
      m.primer_right,
    ];

    const map = new Map((m.blast ?? []).map((b) => [b.db, b]));
    dbs.forEach((db) => {
      const b = map.get(db);
      row.push(b?.amplicon_count ?? null);
      row.push(b?.quality ?? "");
      row.push(b?.top_subject ?? "");
      row.push(b?.top_start ?? null);
      row.push(b?.top_end ?? null);
      row.push(b?.gene_label ?? "");
    });

    return row;
  });

  const summary: Array<Array<string | number>> = [
    ["created", new Date().toLocaleString()],
    ["ref_db", res.ref_db],
    ["ref_entry", res.ref_entry],
    ["ref_start", res.ref_start],
    ["ref_end", res.ref_end],
    ["ref_length", res.ref_length],
    ["alt_db", res.alt_db],
    ["alt_entry", res.alt_entry],
    ["alt_start", res.alt_start],
    ["alt_end", res.alt_end],
    ["alt_strand", res.alt_strand],
    ["alt_length", res.alt_length],
    ["mapped_by_blast", res.mapped_by_blast ? "true" : "false"],
    ["primer_pairs_generated", res.primer_pairs_generated],
    ["markers", res.markers.length],
  ];

  return [
    { name: "Markers", data: [header, ...rows] },
    { name: "Summary", data: summary },
  ];
};

export const CapsPrimerPanel: React.FC = () => {
  const { showToast } = useToast();
  const { options: localDbOptions } = useLocalBlastDbOptions();
  const localDbOptionsWithCustom = useMemo(
    () => withCustomDbOption(localDbOptions, "手動入力（任意の makeblastdb パス）"),
    [localDbOptions],
  );

  const [refDbChoice, setRefDbChoice] = useState<string>(() => {
    const stored = loadPreferredLocalDbPaths() ?? [];
    const real = stored.filter((p) => p && p !== query_TO_ref_VIRTUAL_DB_VALUE);
    return real[0] ?? FALLBACK_LOCAL_DB_OPTIONS[0]?.value ?? "";
  });
  const [refDbCustom, setRefDbCustom] = useState<string>("");
  const refDb = useMemo(
    () => (refDbChoice === CUSTOM_DB_VALUE ? normalizeLocalDbValue(refDbCustom) : refDbChoice),
    [refDbChoice, refDbCustom],
  );

  const [altDbChoice, setAltDbChoice] = useState<string>(() => {
    const stored = loadPreferredLocalDbPaths() ?? [];
    const real = stored.filter((p) => p && p !== query_TO_ref_VIRTUAL_DB_VALUE);
    return real[1] ?? real[0] ?? FALLBACK_LOCAL_DB_OPTIONS[1]?.value ?? FALLBACK_LOCAL_DB_OPTIONS[0]?.value ?? "";
  });
  const [altDbCustom, setAltDbCustom] = useState<string>("");
  const altDb = useMemo(
    () => (altDbChoice === CUSTOM_DB_VALUE ? normalizeLocalDbValue(altDbCustom) : altDbChoice),
    [altDbChoice, altDbCustom],
  );

  useEffect(() => {
    if (!localDbOptionsWithCustom.length) return;
    const valid = new Set(localDbOptionsWithCustom.map((o) => o.value));
    if (!refDbChoice || !valid.has(refDbChoice)) {
      setRefDbChoice(localDbOptionsWithCustom[0]?.value ?? "");
    }
    if (!altDbChoice || !valid.has(altDbChoice)) {
      setAltDbChoice(localDbOptionsWithCustom[Math.min(1, localDbOptionsWithCustom.length - 1)]?.value ?? "");
    }
  }, [altDbChoice, localDbOptionsWithCustom, refDbChoice]);

  useEffect(() => {
    const hasref = localDbOptionsWithCustom.some((o) => o.value === "UserDB_ref");
    const hasALT = localDbOptionsWithCustom.some((o) => o.value === "UserDB_ALT");
    if (hasref) setChrMapRefDb((prev) => prev || "UserDB_ref");
    else if (hasALT) setChrMapRefDb((prev) => prev || "UserDB_ALT");
    else if (localDbOptionsWithCustom[0]?.value) setChrMapRefDb((prev) => prev || localDbOptionsWithCustom[0].value);
  }, [localDbOptionsWithCustom]);

  const [refEntry, setRefEntry] = useState<string>("");
  const [refStart, setRefStart] = useState<number>(1);
  const [refEnd, setRefEnd] = useState<number>(5000);
  const [refChroms, setRefChroms] = useState<BlastDbChromosome[] | null>(null);

  const [mapAltByBlast, setMapAltByBlast] = useState<boolean>(true);
  const [altEntry, setAltEntry] = useState<string>("");
  const [altStart, setAltStart] = useState<number>(1);
  const [altEnd, setAltEnd] = useState<number>(5000);
  const [altStrand, setAltStrand] = useState<"plus" | "minus">("plus");
  const [altChroms, setAltChroms] = useState<BlastDbChromosome[] | null>(null);

  const [chrMapJobId, setChrMapJobId] = useState<string | null>(null);
  const [chrMapJobInfo, setChrMapJobInfo] = useState<JobInfo | null>(null);
  const [chrMapRefDb, setChrMapRefDb] = useState<string>("UserDB_ref");

  const [productMin, setProductMin] = useState<number>(200);
  const [productMax, setProductMax] = useState<number>(800);
  const [primerNumReturn, setPrimerNumReturn] = useState<number>(200);
  const [maxMarkers, setMaxMarkers] = useState<number>(200);

  const [enzymesText, setEnzymesText] = useState<string>(DEFAULT_ENZYMES_TEXT);
  const enzymes = useMemo(() => normalizeEnzymeList(enzymesText), [enzymesText]);
  const [enzymesPerPrimer, setEnzymesPerPrimer] = useState<number>(2);
  const [maxCutsPerAllele, setMaxCutsPerAllele] = useState<number>(3);
  const [minFragmentLen, setMinFragmentLen] = useState<number>(30);
  const [requirePerfectPrimersInAlt, setRequirePerfectPrimersInAlt] = useState<boolean>(true);

  const [blastCheckDbs, setBlastCheckDbs] = useState<string[]>(() => {
    const stored = loadPreferredLocalDbPaths() ?? [];
    const real = stored.filter((p) => p && p !== query_TO_ref_VIRTUAL_DB_VALUE);
    return real.length ? real : FALLBACK_LOCAL_DB_OPTIONS.map((o) => o.value);
  });
  const [blastCheckCustomDb, setBlastCheckCustomDb] = useState<string>("");
  const effectiveBlastCheckDbs = useMemo(() => {
    const manual = normalizeLocalDbValue(blastCheckCustomDb);
    const list = [...blastCheckDbs, ...(manual ? [manual] : [])];
    return Array.from(new Set(list)).filter(Boolean);
  }, [blastCheckCustomDb, blastCheckDbs]);
  const [blastMaxTargetSeqs, setBlastMaxTargetSeqs] = useState<number>(25);
  const [blastNumThreads, setBlastNumThreads] = useState<number | null>(null);

  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<CapsDesignResponse | null>(null);
  const [jobId, setJobId] = useState<string | null>(null);
  const [jobInfo, setJobInfo] = useState<JobInfo | null>(null);

  const STORAGE_KEY = "seq_workbench_caps_panel";

  useEffect(() => {
    if (typeof window === "undefined") return;
    try {
      const raw = window.localStorage.getItem(STORAGE_KEY);
      if (!raw) return;
      const saved = JSON.parse(raw) as Partial<{
        refDbChoice: string;
        refDbCustom: string;
        altDbChoice: string;
        altDbCustom: string;
        refEntry: string;
        refStart: number;
        refEnd: number;
        mapAltByBlast: boolean;
        altEntry: string;
        altStart: number;
        altEnd: number;
        altStrand: "plus" | "minus";
        productMin: number;
        productMax: number;
        primerNumReturn: number;
        maxMarkers: number;
        enzymesText: string;
        enzymesPerPrimer: number;
        maxCutsPerAllele: number;
        minFragmentLen: number;
        requirePerfectPrimersInAlt: boolean;
        blastCheckDbs: string[];
        blastCheckCustomDb: string;
        blastMaxTargetSeqs: number;
        blastNumThreads: number | null;
      }>;
      if (typeof saved.refDbChoice === "string") setRefDbChoice(normalizeLocalDbValue(saved.refDbChoice));
      if (typeof saved.refDbCustom === "string") setRefDbCustom(saved.refDbCustom);
      if (typeof saved.altDbChoice === "string") setAltDbChoice(normalizeLocalDbValue(saved.altDbChoice));
      if (typeof saved.altDbCustom === "string") setAltDbCustom(saved.altDbCustom);
      if (typeof saved.refEntry === "string") setRefEntry(saved.refEntry);
      if (typeof saved.refStart === "number") setRefStart(Math.max(1, Math.floor(saved.refStart)));
      if (typeof saved.refEnd === "number") setRefEnd(Math.max(1, Math.floor(saved.refEnd)));
      if (typeof saved.mapAltByBlast === "boolean") setMapAltByBlast(saved.mapAltByBlast);
      if (typeof saved.altEntry === "string") setAltEntry(saved.altEntry);
      if (typeof saved.altStart === "number") setAltStart(Math.max(1, Math.floor(saved.altStart)));
      if (typeof saved.altEnd === "number") setAltEnd(Math.max(1, Math.floor(saved.altEnd)));
      if (saved.altStrand === "plus" || saved.altStrand === "minus") setAltStrand(saved.altStrand);
      if (typeof saved.productMin === "number") setProductMin(Math.max(50, Math.floor(saved.productMin)));
      if (typeof saved.productMax === "number") setProductMax(Math.max(50, Math.floor(saved.productMax)));
      if (typeof saved.primerNumReturn === "number") setPrimerNumReturn(Math.max(1, Math.floor(saved.primerNumReturn)));
      if (typeof saved.maxMarkers === "number") setMaxMarkers(Math.max(1, Math.floor(saved.maxMarkers)));
      if (typeof saved.enzymesText === "string") setEnzymesText(saved.enzymesText);
      if (typeof saved.enzymesPerPrimer === "number") setEnzymesPerPrimer(Math.max(1, Math.floor(saved.enzymesPerPrimer)));
      if (typeof saved.maxCutsPerAllele === "number") setMaxCutsPerAllele(Math.max(0, Math.floor(saved.maxCutsPerAllele)));
      if (typeof saved.minFragmentLen === "number") setMinFragmentLen(Math.max(1, Math.floor(saved.minFragmentLen)));
      if (typeof saved.requirePerfectPrimersInAlt === "boolean") setRequirePerfectPrimersInAlt(saved.requirePerfectPrimersInAlt);
      if (Array.isArray(saved.blastCheckDbs) && saved.blastCheckDbs.length) {
        const normalized = Array.from(
          new Set(
            saved.blastCheckDbs
              .filter((v) => typeof v === "string")
              .map((v) => normalizeLocalDbValue(v))
              .filter(Boolean),
          ),
        );
        if (normalized.length) setBlastCheckDbs(normalized);
      }
      if (typeof saved.blastCheckCustomDb === "string") setBlastCheckCustomDb(saved.blastCheckCustomDb);
      if (typeof saved.blastMaxTargetSeqs === "number") setBlastMaxTargetSeqs(Math.max(1, Math.floor(saved.blastMaxTargetSeqs)));
      if (typeof saved.blastNumThreads === "number" || saved.blastNumThreads === null) setBlastNumThreads(saved.blastNumThreads ?? null);
    } catch {
      // ignore
    }
  }, []);

  useEffect(() => {
    if (typeof window === "undefined") return;
    const payload = {
      refDbChoice,
      refDbCustom,
      altDbChoice,
      altDbCustom,
      refEntry,
      refStart,
      refEnd,
      mapAltByBlast,
      altEntry,
      altStart,
      altEnd,
      altStrand,
      productMin,
      productMax,
      primerNumReturn,
      maxMarkers,
      enzymesText,
      enzymesPerPrimer,
      maxCutsPerAllele,
      minFragmentLen,
      requirePerfectPrimersInAlt,
      blastCheckDbs,
      blastCheckCustomDb,
      blastMaxTargetSeqs,
      blastNumThreads,
    };
    try {
      window.localStorage.setItem(STORAGE_KEY, JSON.stringify(payload));
    } catch {
      // ignore
    }
  }, [
    altDbChoice,
    altDbCustom,
    altEnd,
    altEntry,
    altStart,
    altStrand,
    blastCheckCustomDb,
    blastCheckDbs,
    blastMaxTargetSeqs,
    blastNumThreads,
    enzymesPerPrimer,
    enzymesText,
    mapAltByBlast,
    maxCutsPerAllele,
    maxMarkers,
    minFragmentLen,
    primerNumReturn,
    productMax,
    productMin,
    refDbChoice,
    refDbCustom,
    refEnd,
    refEntry,
    refStart,
    requirePerfectPrimersInAlt,
  ]);

  useEffect(() => {
    let cancelled = false;
    const db = refDb.trim();
    if (!db) {
      setRefChroms(null);
      return;
    }
    bioapiClient
      .listDbChromosomes(db)
      .then((rows) => {
        if (cancelled) return;
        setRefChroms(rows);
      })
      .catch(() => {
        if (cancelled) return;
        setRefChroms(null);
      });
    return () => {
      cancelled = true;
    };
  }, [refDb]);

  useEffect(() => {
    let cancelled = false;
    const db = altDb.trim();
    if (!db) {
      setAltChroms(null);
      return;
    }
    bioapiClient
      .listDbChromosomes(db)
      .then((rows) => {
        if (cancelled) return;
        setAltChroms(rows);
      })
      .catch(() => {
        if (cancelled) return;
        setAltChroms(null);
      });
    return () => {
      cancelled = true;
    };
  }, [altDb]);

  useEffect(() => {
    if (!jobId) return;
    const ctrl = new AbortController();

    const runPoll = async () => {
      try {
        const info = await pollJobUntilDone(jobId, {
          onUpdate: (i) => setJobInfo(i),
          intervalMs: 900,
          signal: ctrl.signal,
        });
        if (ctrl.signal.aborted) return;
        if (info.status !== "succeeded") {
          const msg = info.error ?? "CAPS 生成に失敗しました。";
          setError(msg);
          showToast(msg, info.status === "canceled" ? "info" : "error");
          return;
        }
        const res = await bioapiClient.getJobResult<CapsDesignResponse>(jobId);
        if (ctrl.signal.aborted) return;
        setResult(res);
        showToast("CAPS 生成が完了しました", "success");
      } catch (e) {
        if (ctrl.signal.aborted) return;
        const msg = e instanceof Error ? e.message : "ジョブ取得に失敗しました。";
        setError(msg);
        showToast(msg, "error");
      } finally {
        if (!ctrl.signal.aborted) {
          setLoading(false);
          setJobId(null);
          setJobInfo(null);
        }
      }
    };

    void runPoll();
    return () => ctrl.abort();
  }, [jobId, showToast]);

  const toggleBlastCheckDb = (path: string) => {
    setBlastCheckDbs((prev) => (prev.includes(path) ? prev.filter((p) => p !== path) : [...prev, path]));
  };

  const run = async () => {
    if (loading || jobId) return;
    if (!refDb) {
      setError("参照 DB を指定してください。");
      return;
    }
    if (!altDb) {
      setError("比較 DB を指定してください。");
      return;
    }
    if (!refEntry.trim()) {
      setError("参照 entry（染色体/コンティグ）を入力してください。");
      return;
    }
    if (!mapAltByBlast && !altEntry.trim()) {
      setError("比較 entry（手動指定）を入力してください。");
      return;
    }
    if (refStart < 1 || refEnd < 1) {
      setError("座標は 1 以上を指定してください。");
      return;
    }
    if (productMin > productMax) {
      setError("product_min <= product_max になるようにしてください。");
      return;
    }
    if (!enzymes.length) {
      setError("制限酵素が空です。デフォルトを入れるか、酵素名を入力してください。");
      return;
    }

    setLoading(true);
    setError(null);
    setResult(null);
    setJobInfo(null);
    try {
      const job = await bioapiClient.createCapsDesignJob({
        ref_db: refDb,
        // UI は chr1 等で選べるが、バックエンドへは実体 seqid を渡す（DB により 1LG6 等に変換）
        ref_entry: resolvedRefEntry,
        ref_start: Math.min(refStart, refEnd),
        ref_end: Math.max(refStart, refEnd),
        alt_db: altDb,
        map_alt_by_blast: mapAltByBlast,
        alt_entry: mapAltByBlast ? null : resolvedAltEntry,
        alt_start: mapAltByBlast ? null : Math.min(altStart, altEnd),
        alt_end: mapAltByBlast ? null : Math.max(altStart, altEnd),
        alt_strand: mapAltByBlast ? "plus" : altStrand,
        product_min: productMin,
        product_max: productMax,
        primer_num_return: primerNumReturn,
        max_markers: maxMarkers,
        enzymes,
        enzymes_per_primer: enzymesPerPrimer,
        max_cuts_per_allele: maxCutsPerAllele,
        min_fragment_len: minFragmentLen,
        require_perfect_primers_in_alt: requirePerfectPrimersInAlt,
        blast_check_dbs: effectiveBlastCheckDbs,
        blast_max_target_seqs: blastMaxTargetSeqs,
        blast_num_threads: blastNumThreads,
      });
      setJobId(job.job_id);
      showToast("CAPS ジョブを開始しました", "info");
    } catch (e) {
      const msg = e instanceof Error ? e.message : "CAPS 生成中にエラーが発生しました。";
      setError(msg);
      showToast(msg, "error");
      setLoading(false);
    } finally {
      // loading は job の完了まで継続
    }
  };

  const cancel = async () => {
    if (!jobId) return;
    try {
      await bioapiClient.cancelJob(jobId);
      showToast("キャンセルを要求しました（実行中の処理はすぐ止まらないことがあります）", "info");
    } catch (e) {
      const msg = e instanceof Error ? e.message : "キャンセルに失敗しました。";
      showToast(msg, "error");
    }
  };

  const buildChromAliasesFor = async (targetDb: string, which: "ref" | "alt") => {
    const db = (targetDb || "").trim();
    const ref = (chrMapRefDb || "").trim();
    if (!db) return;
    if (chrMapJobId) return;
    setChrMapJobId(null);
    setChrMapJobInfo(null);
    try {
      const job = await bioapiClient.createBuildChromAliasesJob({
        db,
        ref_db: ref || "UserDB_ref",
      });
      setChrMapJobId(job.job_id);
      const info = await pollJobUntilDone(job.job_id, {
        onUpdate: (i) => setChrMapJobInfo(i),
        intervalMs: 900,
      });
      if (info.status !== "succeeded") {
        throw new Error(info.error ?? "染色体推定に失敗しました。");
      }
      const rows = await bioapiClient.listDbChromosomes(db);
      if (which === "ref") setRefChroms(rows);
      else setAltChroms(rows);
      showToast("染色体（chr→entry）推定が完了しました", "success");
    } catch (e) {
      const msg = e instanceof Error ? e.message : "染色体推定に失敗しました。";
      showToast(msg, "error");
    } finally {
      setChrMapJobId(null);
      setChrMapJobInfo(null);
    }
  };

  const refChromValue = useMemo(() => {
    const n = normalizeChr(refEntry);
    if (!n) return "";
    const chromSet = new Set((refChroms ?? []).map((c) => c.chrom));
    const fallbackSet = new Set(["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"]);
    return chromSet.has(n) || fallbackSet.has(n) ? n : "";
  }, [refChroms, refEntry]);

  const altChromValue = useMemo(() => {
    const n = normalizeChr(altEntry);
    if (!n) return "";
    const chromSet = new Set((altChroms ?? []).map((c) => c.chrom));
    const fallbackSet = new Set(["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"]);
    return chromSet.has(n) || fallbackSet.has(n) ? n : "";
  }, [altChroms, altEntry]);

  const resolvedRefEntry = useMemo(() => resolveEntryFromChrom(refEntry, refChroms), [refChroms, refEntry]);
  const resolvedAltEntry = useMemo(() => resolveEntryFromChrom(altEntry, altChroms), [altChroms, altEntry]);

  const tsv = useMemo(() => (result ? buildTsv(result) : ""), [result]);
  const md = useMemo(() => (result ? buildMarkdown(result) : ""), [result]);

  const copyText = async (text: string, label: string) => {
    const t = text.trim();
    if (!t) return;
    try {
      await navigator.clipboard.writeText(t);
      showToast(`${label} をコピーしました`, "success");
    } catch {
      showToast(`${label} のコピーに失敗しました`, "error");
    }
  };

  const summary = useMemo(() => {
    if (!result) return null;
    return {
      ref: `${labelForDb(result.ref_db)} / ${result.ref_entry}:${result.ref_start.toLocaleString()}–${result.ref_end.toLocaleString()} (${result.ref_length.toLocaleString()} bp)`,
      alt: `${labelForDb(result.alt_db)} / ${result.alt_entry}:${result.alt_start.toLocaleString()}–${result.alt_end.toLocaleString()} (${result.alt_length.toLocaleString()} bp) ${result.alt_strand} / mapped_by_blast=${result.mapped_by_blast ? "true" : "false"}`,
    };
  }, [result]);

  const renderRow = (m: CapsMarkerRow) => {
    const url = ensemblGeneUrl(m.gene_label ?? null);
    return (
      <tr key={`${m.index}-${m.enzyme}-${m.primer_left}-${m.primer_right}`}>
        <td>{m.index}</td>
        <td>{m.enzyme}</td>
        <td className="blast-desc" style={{ whiteSpace: "normal" }}>
          {m.gene_label ? (
            <EnsemblLinksInline
              geneId={m.gene_label}
              /* We don't have the specific DB for this gene label easily here, but usually it matches the ref/alt DBs or blast DBs.
                 If we omit dbLabel, it infers from ID.
                 If ID is GENE... it goes to Ensembl.
                 If ID is non-standard, it might default to Search.
                 To support Navigator for non-standard IDs, EnsemblLinksInline needs to know if it should use Navigator.
                 The current logic in EnsemblLinksInline:
                 if (localOnly) -> Navigator
                 if (species) -> Ensembl
                 else -> Ensembl Search
                 
                 We verified isLocalOnlyDb checks the DB Label.
                 Here we can pass `res.ref_db` or similar if we assume the gene is from there?
                 But `m.gene_label` might be composite.
                 However, `CapsPrimerPanel` is often used with local DBs.
                 Let's try to pass `refDb` as a hint if it helps, but `m.gene_label` is often just the ID.
                 
                 Actually, looking at `renderRow`:
                 It renders `m.gene_label`.
                 If I use `<EnsemblLinksInline geneId={m.gene_label} dbLabel={refDb} />`, it will respect the Ref DB setting.
                 If Ref DB is "sample_v1" (local only), it will link to Navigator.
                 If Ref DB is "UserDB_ref", it will link to Ensembl.
                 This seems correct for the primary gene label which usually comes from the Ref DB context.
              */
              dbLabel={refDb}
              showLocation={false}
              showTranscript={false}
              showExport={false}
            />
          ) : (
            "-"
          )}
        </td>
        <td>{formatRange(m.ref_product_start, m.ref_product_end)}</td>
        <td>
          {formatRange(m.alt_product_start, m.alt_product_end)} ({m.alt_strand})
        </td>
        <td>
          {m.product_len_ref} / {m.product_len_alt}
        </td>
        <td>
          {formatFragments(m.fragments_ref)} → {formatFragments(m.fragments_alt)}
        </td>
        <td>
          <GelMini
            fragmentsRef={m.fragments_ref}
            fragmentsAlt={m.fragments_alt}
            productLenRef={m.product_len_ref}
            productLenAlt={m.product_len_alt}
          />
        </td>
        <td>{m.mismatch_count}</td>
        <td className="blast-desc">
          <code>{m.primer_left}</code>
          <button type="button" className="link-button" onClick={() => copyText(m.primer_left, "Primer F")}>
            copy
          </button>
        </td>
        <td className="blast-desc">
          <code>{m.primer_right}</code>
          <button type="button" className="link-button" onClick={() => copyText(m.primer_right, "Primer R")}>
            copy
          </button>
        </td>
        <td className="blast-desc" style={{ whiteSpace: "normal" }}>
          {(m.blast ?? []).map((b) => (
            <div key={b.db}>
              <strong>{b.db}</strong>: {b.amplicon_count} {b.quality ? `(${b.quality})` : ""}
              {b.top_subject && b.top_start != null && b.top_end != null ? (
                (() => {
                  return (
                    <>
                      {" "}
                      / {b.top_subject}:{formatRange(b.top_start, b.top_end)}
                      {" "}
                      <EnsemblLinksInline
                        geneId={b.gene_label ?? null}
                        dbLabel={b.db}
                        chrom={b.top_subject ?? null}
                        start={b.top_start ?? null}
                        end={b.top_end ?? null}
                        showGene={false}
                      />
                    </>
                  );
                })()
              ) : null}
              {b.gene_label ? (
                <>
                  {" "}
                  /{" "}
                  {(() => {
                    const localOnly = isLocalOnlyDb(b.db);
                    const url = localOnly
                      ? navigatorGeneUrl({ geneId: b.gene_label, dbLabel: b.db })
                      : ensemblGeneUrl(b.gene_label);
                    return url ? (
                      <a href={url} target="_blank" rel="noreferrer">
                        {b.gene_label}
                      </a>
                    ) : (
                      b.gene_label
                    );
                  })()}
                </>
              ) : null}
            </div>
          ))}
        </td>
      </tr>
    );
  };

  return (
    <section className="seq-result-block">
      <h2 className="panel-title">CAPSプライマー作成（大量生成）</h2>
      <p className="panel-hint">
        指定したゲノム範囲で Primer3 による多数プライマーペアを作り、参照/比較の制限酵素切断パターン差から
        共優勢（CAPS）候補を抽出します。比較側は既定で BLAST で対応領域を推定します。
      </p>

      <div className="primer-row" style={{ marginBottom: "0.4rem" }}>
        <button
          type="button"
          className="seq-button secondary"
          onClick={() => {
            if (!result) return;
            downloadXlsx(buildXlsxSheets(result), "caps_markers");
          }}
          disabled={!result}
        >
          結果を Excel (.xlsx) として保存
        </button>
        <button
          type="button"
          className="seq-button secondary"
          onClick={() => {
            if (!md) return;
            downloadMarkdown(md, "caps_markers");
          }}
          disabled={!result}
        >
          結果を Markdown として保存
        </button>
        <button
          type="button"
          className="seq-button secondary"
          onClick={() => {
            if (!md) return;
            openPrintViewForMarkdown(md, "CAPS プライマーレポート");
          }}
          disabled={!result}
        >
          印刷用ビューを開く（PDF 保存に利用）
        </button>
        <button type="button" className="seq-button secondary" onClick={() => copyText(tsv, "TSV")} disabled={!tsv}>
          TSV をコピー
        </button>
      </div>

      <div className="primer-grid">
        <div className="primer-controls">
          <h3>入力</h3>
          <div className="form-grid">
            <label className="seq-label">
              参照DB:
              <select className="seq-input" value={refDbChoice} onChange={(e) => setRefDbChoice(e.target.value)} disabled={loading}>
                {localDbOptionsWithCustom.map((o) => (
                  <option key={o.value} value={o.value}>
                    {o.label}
                  </option>
                ))}
              </select>
            </label>
            <label className="seq-label">
              参照染色体（簡易）:
              <select className="seq-input" value={refChromValue} onChange={(e) => setRefEntry(e.target.value)} disabled={loading}>
                <option value="">(手入力)</option>
                {(refChroms ?? []).map((c) => (
                  <option key={`${c.chrom}-${c.entry}`} value={c.chrom}>
                    {c.chrom} → {c.entry}
                  </option>
                ))}
                {!refChroms?.length && ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"].map((c) => (
                  <option key={`fallback-${c}`} value={c}>
                    {c}
                  </option>
                ))}
              </select>
            </label>
            {!refChroms?.length ? (
              <div className="seq-hint grid-span-2">
                <div style={{ marginBottom: "0.25rem" }}>
                  このDBは chr→entry が自動推定できません。ref/ALT を参照して best-effort で chr→entry を推定できます（結果はキャッシュされます）。
                </div>
                <div className="primer-row" style={{ alignItems: "center", gap: "0.6rem", flexWrap: "wrap" }}>
                  <label className="seq-hint" style={{ display: "flex", alignItems: "center", gap: "0.4rem" }}>
                    参照DB:
                    <select
                      className="seq-input"
                      value={chrMapRefDb}
                      onChange={(e) => setChrMapRefDb(e.target.value)}
                      disabled={loading}
                      style={{ width: "auto", minWidth: "12rem" }}
                    >
                      {localDbOptionsWithCustom
                        .filter((o) => o.value !== CUSTOM_DB_VALUE)
                        .map((o) => (
                          <option key={o.value} value={o.value}>
                            {o.label}
                          </option>
                        ))}
                    </select>
                  </label>
                  <button
                    type="button"
                    className="seq-button secondary"
                    onClick={() => void buildChromAliasesFor(refDb, "ref")}
                    disabled={loading || !refDb}
                  >
                    参照DBの染色体（chr→entry）を推定
                  </button>
                </div>
                <JobProgressCard title="染色体推定" jobId={chrMapJobId} job={chrMapJobInfo} />
              </div>
            ) : null}
            {refDbChoice === CUSTOM_DB_VALUE && (
              <label className="seq-label grid-span-2">
                参照DBパス:
                <input className="seq-input" type="text" value={refDbCustom} onChange={(e) => setRefDbCustom(e.target.value)} placeholder="/path/to/db_prefix" />
              </label>
            )}
            <label className="seq-label grid-span-2">
              参照 entry（染色体/コンティグ）:
              <input className="seq-input" type="text" value={refEntry} onChange={(e) => setRefEntry(e.target.value)} placeholder="例: chr1（DBに応じて自動変換） / 1LG6 / NC_..." />
            </label>
            {refEntry.trim() && resolvedRefEntry && resolvedRefEntry !== refEntry.trim() && (
              <p className="seq-hint grid-span-2">
                entry 解決: <code>{refEntry.trim()}</code> → <code>{resolvedRefEntry}</code>
              </p>
            )}
            <label className="seq-label">
              start:
              <input className="seq-input" type="number" min={1} value={refStart} onChange={(e) => setRefStart(Math.max(1, Number(e.target.value) || 1))} />
            </label>
            <label className="seq-label">
              end:
              <input className="seq-input" type="number" min={1} value={refEnd} onChange={(e) => setRefEnd(Math.max(1, Number(e.target.value) || 1))} />
            </label>
          </div>

          <hr />

          <div className="form-grid">
            <label className="seq-label">
              比較DB:
              <select className="seq-input" value={altDbChoice} onChange={(e) => setAltDbChoice(e.target.value)} disabled={loading}>
                {localDbOptionsWithCustom.map((o) => (
                  <option key={o.value} value={o.value}>
                    {o.label}
                  </option>
                ))}
              </select>
            </label>
            <label className="seq-hint grid-span-2" style={{ display: "flex", alignItems: "center", gap: "0.4rem" }}>
              <input type="checkbox" checked={mapAltByBlast} onChange={(e) => setMapAltByBlast(e.target.checked)} />
              比較側の対応領域を BLAST で推定する（推奨）
            </label>
            {altDbChoice === CUSTOM_DB_VALUE && (
              <label className="seq-label grid-span-2">
                比較DBパス:
                <input className="seq-input" type="text" value={altDbCustom} onChange={(e) => setAltDbCustom(e.target.value)} placeholder="/path/to/db_prefix" />
              </label>
            )}
            {!mapAltByBlast && (
              <>
                <label className="seq-label">
                  比較染色体（簡易）:
                  <select className="seq-input" value={altChromValue} onChange={(e) => setAltEntry(e.target.value)} disabled={loading}>
                    <option value="">(手入力)</option>
                    {(altChroms ?? []).map((c) => (
                      <option key={`alt-${c.chrom}-${c.entry}`} value={c.chrom}>
                        {c.chrom} → {c.entry}
                      </option>
                    ))}
                    {!altChroms?.length && ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7"].map((c) => (
                      <option key={`alt-fallback-${c}`} value={c}>
                        {c}
                      </option>
                    ))}
                  </select>
                </label>
                {!altChroms?.length ? (
                  <div className="seq-hint grid-span-2">
                    <div style={{ marginBottom: "0.25rem" }}>
                      このDBは chr→entry が自動推定できません。ref/ALT を参照して best-effort で chr→entry を推定できます（結果はキャッシュされます）。
                    </div>
                    <div className="primer-row" style={{ alignItems: "center", gap: "0.6rem", flexWrap: "wrap" }}>
                      <label className="seq-hint" style={{ display: "flex", alignItems: "center", gap: "0.4rem" }}>
                        参照DB:
                        <select
                          className="seq-input"
                          value={chrMapRefDb}
                          onChange={(e) => setChrMapRefDb(e.target.value)}
                          disabled={loading}
                          style={{ width: "auto", minWidth: "12rem" }}
                        >
                          {localDbOptionsWithCustom
                            .filter((o) => o.value !== CUSTOM_DB_VALUE)
                            .map((o) => (
                              <option key={o.value} value={o.value}>
                                {o.label}
                              </option>
                            ))}
                        </select>
                      </label>
                      <button
                        type="button"
                        className="seq-button secondary"
                        onClick={() => void buildChromAliasesFor(altDb, "alt")}
                        disabled={loading || !altDb}
                      >
                        比較DBの染色体（chr→entry）を推定
                      </button>
                    </div>
                    <JobProgressCard title="染色体推定" jobId={chrMapJobId} job={chrMapJobInfo} />
                  </div>
                ) : null}
                <label className="seq-label grid-span-2">
                  比較 entry（手動）:
                  <input className="seq-input" type="text" value={altEntry} onChange={(e) => setAltEntry(e.target.value)} />
                </label>
                {altEntry.trim() && resolvedAltEntry && resolvedAltEntry !== altEntry.trim() && (
                  <p className="seq-hint grid-span-2">
                    entry 解決: <code>{altEntry.trim()}</code> → <code>{resolvedAltEntry}</code>
                  </p>
                )}
                <label className="seq-label">
                  alt start:
                  <input className="seq-input" type="number" min={1} value={altStart} onChange={(e) => setAltStart(Math.max(1, Number(e.target.value) || 1))} />
                </label>
                <label className="seq-label">
                  alt end:
                  <input className="seq-input" type="number" min={1} value={altEnd} onChange={(e) => setAltEnd(Math.max(1, Number(e.target.value) || 1))} />
                </label>
                <label className="seq-label">
                  alt strand:
                  <select className="seq-input" value={altStrand} onChange={(e) => setAltStrand(e.target.value as "plus" | "minus")}>
                    <option value="plus">plus</option>
                    <option value="minus">minus</option>
                  </select>
                </label>
              </>
            )}
          </div>

          <hr />

          <div className="form-grid">
            <label className="seq-label">
              product_min:
              <input className="seq-input" type="number" min={50} value={productMin} onChange={(e) => setProductMin(Math.max(50, Number(e.target.value) || 200))} />
            </label>
            <label className="seq-label">
              product_max:
              <input className="seq-input" type="number" min={50} value={productMax} onChange={(e) => setProductMax(Math.max(50, Number(e.target.value) || 800))} />
            </label>
            <label className="seq-label">
              primer_num_return:
              <input className="seq-input" type="number" min={1} max={2000} value={primerNumReturn} onChange={(e) => setPrimerNumReturn(Math.max(1, Math.min(2000, Number(e.target.value) || 200)))} />
            </label>
            <label className="seq-label">
              max_markers:
              <input className="seq-input" type="number" min={1} max={5000} value={maxMarkers} onChange={(e) => setMaxMarkers(Math.max(1, Math.min(5000, Number(e.target.value) || 200)))} />
            </label>
            <label className="seq-hint grid-span-2" style={{ display: "flex", alignItems: "center", gap: "0.4rem" }}>
              <input type="checkbox" checked={requirePerfectPrimersInAlt} onChange={(e) => setRequirePerfectPrimersInAlt(e.target.checked)} />
              比較側でもプライマー結合部が完全一致する候補のみ返す（安全）
            </label>
          </div>

          <div className="form-grid" style={{ marginTop: "0.35rem" }}>
            <label className="seq-label grid-span-2">
              制限酵素（Biopython名, 改行/スペース区切り）:
              <textarea className="seq-textarea" value={enzymesText} onChange={(e) => setEnzymesText(e.target.value)} rows={10} />
              <div className="primer-row" style={{ marginTop: "0.25rem" }}>
                <button type="button" className="seq-button secondary" onClick={() => setEnzymesText(DEFAULT_ENZYMES_TEXT)} disabled={loading}>
                  デフォルトを入れる
                </button>
                <span className="seq-hint">現在 {enzymes.length} 個</span>
              </div>
            </label>
            <label className="seq-label">
              enzymes_per_primer:
              <input className="seq-input" type="number" min={1} max={20} value={enzymesPerPrimer} onChange={(e) => setEnzymesPerPrimer(Math.max(1, Math.min(20, Number(e.target.value) || 2)))} />
            </label>
            <label className="seq-label">
              max_cuts:
              <input className="seq-input" type="number" min={0} max={20} value={maxCutsPerAllele} onChange={(e) => setMaxCutsPerAllele(Math.max(0, Math.min(20, Number(e.target.value) || 3)))} />
            </label>
            <label className="seq-label">
              min_fragment_len:
              <input className="seq-input" type="number" min={1} max={500} value={minFragmentLen} onChange={(e) => setMinFragmentLen(Math.max(1, Math.min(500, Number(e.target.value) || 30)))} />
            </label>
          </div>

          <hr />
          <h4 style={{ margin: "0.4rem 0 0.2rem" }}>一意性チェック（ローカル BLAST）</h4>
          <div className="form-grid">
            <div className="tag-row grid-span-2" style={{ flexWrap: "wrap", gap: "0.6rem", marginBottom: "0.1rem" }}>
              {localDbOptions.map((opt) => (
                <label key={opt.value} className="seq-hint" style={{ display: "flex", alignItems: "center", gap: "0.3rem" }}>
                  <input type="checkbox" checked={blastCheckDbs.includes(opt.value)} onChange={() => toggleBlastCheckDb(opt.value)} />
                  {opt.label}
                </label>
              ))}
            </div>
            <label className="seq-label grid-span-2">
              追加DB（任意）:
              <input className="seq-input" type="text" value={blastCheckCustomDb} onChange={(e) => setBlastCheckCustomDb(e.target.value)} placeholder="/path/to/db_prefix" />
            </label>
            <label className="seq-label">
              max_target_seqs:
              <input className="seq-input" type="number" min={1} max={200} value={blastMaxTargetSeqs} onChange={(e) => setBlastMaxTargetSeqs(Math.max(1, Math.min(200, Number(e.target.value) || 25)))} />
            </label>
            <label className="seq-label">
              num_threads:
              <input
                className="seq-input"
                type="number"
                min={1}
                value={blastNumThreads ?? ""}
                onChange={(e) => {
                  const v = Number(e.target.value);
                  setBlastNumThreads(Number.isNaN(v) ? null : Math.max(1, Math.floor(v)));
                }}
                placeholder="自動"
              />
            </label>
          </div>

          <button type="button" className="seq-button" onClick={run} disabled={loading}>
            {loading ? "実行中..." : "CAPS マーカーを作る"}
          </button>
          <JobProgressCard
            title="CAPS 生成"
            jobId={jobId}
            job={jobInfo}
            onCancel={jobId ? cancel : null}
            cancelDisabled={!jobId}
          />
          {error && <p className="seq-error">エラー: {error}</p>}
        </div>

        <div className="primer-results">
          <h3>結果</h3>
          {!result && !loading && <p className="seq-hint">条件を入力して「CAPS マーカーを作る」を押してください。</p>}
          {summary && (
            <div className="seq-result-block">
              <div className="tag-row">
                <span className="tag-label">ref</span>
                <code className="tag-db">{summary.ref}</code>
              </div>
              <div className="tag-row">
                <span className="tag-label">alt</span>
                <code className="tag-db">{summary.alt}</code>
              </div>
              {result?.warnings?.length ? (
                <p className="seq-hint" style={{ marginTop: "0.4rem" }}>
                  warnings: {result.warnings.join(" / ")}
                </p>
              ) : null}
            </div>
          )}

          {result && (
            <>
              <p className="seq-hint">
                markers: {result.markers.length.toLocaleString()} / primer3 pairs:{" "}
                {result.primer_pairs_generated.toLocaleString()}
              </p>
              {result.markers.length ? (
                <div className="table-scroll">
                  <table className="seq-table">
                    <thead>
                      <tr>
                        <th>#</th>
                        <th>Enzyme</th>
                        <th>Gene</th>
                        <th>Ref range</th>
                        <th>Alt range</th>
                        <th>len(ref/alt)</th>
                        <th>Fragments</th>
                        <th>Gel</th>
                        <th>mismatch</th>
                        <th>Primer F</th>
                        <th>Primer R</th>
                        <th>BLAST check</th>
                      </tr>
                    </thead>
                    <tbody>{result.markers.map(renderRow)}</tbody>
                  </table>
                </div>
              ) : (
                <p className="seq-hint">
                  条件を満たす CAPS 候補が見つかりませんでした。範囲/制限酵素/marker数/完全一致条件を調整してください。
                </p>
              )}
            </>
          )}
        </div>
      </div>
    </section>
  );
};

