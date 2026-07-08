import React, { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { bioapiClient } from "../api/bioapiClient";
import type { PrimerDesignResponse, PrimerPair } from "../types/primers";
import type { BlastResponse, BlastHit, NCBITarget } from "../types/blast";
import type { JobInfo } from "../types/jobs";
import { computePrimerAmplicons, countLocalHits } from "../utils/primerBlast";
import { downloadMarkdown, openPrintViewForMarkdown } from "../utils/exportReport";
import { runBlastBatchLocalJob } from "../utils/blastBatchLocalJob";
import { useWorkbench } from "../utils/workbenchContext";
import { ensemblGeneUrl } from "../utils/ensembl";
import { FeatureSequenceView, type PrimerRange1Based, type Range1Based } from "./FeatureSequenceView";
import { JobProgressCard } from "./JobProgressCard";
import {
  DEFAULT_BLAST_DB_BASE,
  labelForDbPath,
  relabelLocalBlastHits,
  useLocalBlastDbOptions,
  usePreferredLocalDbPaths,
  normalizeLocalDbValue,
  withCustomDbOption,
} from "../utils/localBlastDbs";
import { useToast } from "./ToastProvider";
import { useLocalBlastMode } from "../utils/localBlastMode";

interface ExonPrimerPanelProps {
  mode?: "exon" | "sequence";
  initialSequencingMode?: boolean;
  initialSeqProductMin?: number;
  initialSeqProductMax?: number;
  initialSeqOverlap?: number;
  initialSeqMargin?: number;
}

// CDS/エキソンをハイライトしながら、サブ領域指定でプライマー設計+BLASTを行うタブ
// mode="exon"  : 従来の CDS/エキソン増幅タブ
// mode="sequence": シーケンスプライマー用タブ（SequencePrimerPanel から利用）
export const ExonPrimerPanel: React.FC<ExonPrimerPanelProps> = ({
  mode = "exon",
  initialSequencingMode,
  initialSeqProductMin,
  initialSeqProductMax,
  initialSeqOverlap,
  initialSeqMargin,
}) => {
  const { showToast } = useToast();
  const { options: localDbOptions } = useLocalBlastDbOptions();
  const [loadGeneDb, setLoadGeneDb] = useState<string>("auto");
  const [loadGeneCustomDb, setLoadGeneCustomDb] = useState<string>("");
  const [loadGeneName, setLoadGeneName] = useState<string>("");
  const [loadGeneMargin, setLoadGeneMargin] = useState<number | "">("");
  const [loadGeneLoading, setLoadGeneLoading] = useState<boolean>(false);

  const [sequence, setSequence] = useState<string>("");
  const [exonRanges, setExonRanges] = useState<string>(""); // 例: 100-200,300-350
  const [cdsRange, setCdsRange] = useState<string>(""); // 例: 150-900
  const [geneId, setGeneId] = useState<string>("");
  const [species, setSpecies] = useState<string>("");
  const [structureLoading, setStructureLoading] = useState<boolean>(false);
  const [structureError, setStructureError] = useState<string | null>(null);

  const [subStart, setSubStart] = useState<number | null>(null);
  const [subEnd, setSubEnd] = useState<number | null>(null);

  const [productSizeRange, setProductSizeRange] = useState<string>("100-400");
  const [numReturn, setNumReturn] = useState<number>(5);

  // PCR 条件（主に Tm）。通常は隠しておき、必要なときだけ展開する。
  const [showPcrOptions, setShowPcrOptions] = useState<boolean>(false);
  const [pcrOptTm, setPcrOptTm] = useState<number>(60.0);
  const [pcrMinTm, setPcrMinTm] = useState<number>(57.0);
  const [pcrMaxTm, setPcrMaxTm] = useState<number>(63.0);
  const [pcrMinSize, setPcrMinSize] = useState<number>(18);
  const [pcrOptSize, setPcrOptSize] = useState<number>(20);
  const [pcrMaxSize, setPcrMaxSize] = useState<number>(27);
  const [pcrMinGc, setPcrMinGc] = useState<number>(20.0);
  const [pcrMaxGc, setPcrMaxGc] = useState<number>(80.0);
  const [pcrSaltMonovalent, setPcrSaltMonovalent] = useState<number>(50.0);
  const [pcrDnaConc, setPcrDnaConc] = useState<number>(50.0);

  const [designResult, setDesignResult] = useState<PrimerDesignResponse | null>(null);
  const [selectedPair, setSelectedPair] = useState<number | null>(null);
  const [selectedPairIndices, setSelectedPairIndices] = useState<number[]>([]);
  const [designError, setDesignError] = useState<string | null>(null);
  const [designLoading, setDesignLoading] = useState<boolean>(false);
  // Sequencing Mode States
  const [sequencingMode, setSequencingMode] = useState<boolean>(
    initialSequencingMode ?? false,
  );
  const [seqProductMin, setSeqProductMin] = useState<number>(
    initialSeqProductMin ?? 500,
  );
  const [seqProductMax, setSeqProductMax] = useState<number>(
    initialSeqProductMax ?? 800,
  );
  const [seqOverlap, setSeqOverlap] = useState<number>(
    initialSeqOverlap ?? 50,
  );
  const [seqMargin, setSeqMargin] = useState<number>(
    initialSeqMargin ?? 50,
  );

  const [selectedLocalDbs, setSelectedLocalDbs] = usePreferredLocalDbPaths();
  const [customLocalDb, setCustomLocalDb] = useState<string>("");
  const effectiveLocalDbs = useMemo(() => {
    const manual = normalizeLocalDbValue(customLocalDb);
    const list = [...selectedLocalDbs, ...(manual ? [manual] : [])];
    return Array.from(new Set(list)).filter(Boolean);
  }, [customLocalDb, selectedLocalDbs]);
  const [blastUseLocal, setBlastUseLocal] = useState<boolean>(true);
  const [blastUseNcbi, setBlastUseNcbi] = useState<boolean>(false);
  const [ncbiOrganism, setNcbiOrganism] = useState<boolean>(false);
  const [ncbiArabidopsis, setNcbiArabidopsis] = useState<boolean>(false);
  const [ncbiCustomQuery, setNcbiCustomQuery] = useState<string>("");
  const [blastMaxHits, setBlastMaxHits] = useState<number>(5);
  const [blastTask, setBlastTask] = useState<string>("blastn-short");
  const [blastEvalue, setBlastEvalue] = useState<number>(1e-5);
  const [blastMaxHsps, setBlastMaxHsps] = useState<number | null>(null);
  const [blastNumThreads, setBlastNumThreads] = useState<number | null>(null);
  const [localMode, setLocalMode] = useLocalBlastMode();
  const [blastError, setBlastError] = useState<string | null>(null);
  const [blastLoading, setBlastLoading] = useState<boolean>(false);
  const [blastJobId, setBlastJobId] = useState<string | null>(null);
  const [blastJobInfo, setBlastJobInfo] = useState<JobInfo | null>(null);
  const [blastLeft, setBlastLeft] = useState<BlastResponse | null>(null);
  const [blastRight, setBlastRight] = useState<BlastResponse | null>(null);
  const [primerBlastResults, setPrimerBlastResults] = useState<
    Array<{ pairIndex: number; left: BlastResponse; right: BlastResponse }>
  >([]);
  const [selectionRange, setSelectionRange] = useState<{ start: number; end: number } | null>(
    null,
  );
  const [geneGuess, setGeneGuess] = useState<string>("");
  const [cdsSequenceInput, setCdsSequenceInput] = useState<string>("");
  const [exonSequencesInput, setExonSequencesInput] = useState<string>("");
  const [infoMessage, setInfoMessage] = useState<string>("");
  const [fastaInput, setFastaInput] = useState<string>("");
  const [genomeHeader, setGenomeHeader] = useState<string>("genome_sequence");
  const fastaHighlightRef = useRef<HTMLDivElement | null>(null);
  const [guideRegion, setGuideRegion] = useState<{ start: number; end: number; ts: number } | null>(
    null,
  );
  const [primerPatternTab, setPrimerPatternTab] = useState<"pattern1" | "pattern2">("pattern1");
  const [exonPrimerPlans, setExonPrimerPlans] = useState<
    Array<{
      kind?: "exon" | "cds";
      exonIndex: number;
      exonEndIndex?: number;
      range: [number, number];
      result?: PrimerDesignResponse;
      subIndex?: number;
      error?: string;
      loading?: boolean;
      filteredCandidates?: Array<{
        pair: PrimerPair;
        left?: BlastResponse;
        right?: BlastResponse;
        note?: string;
        ampliconCount?: number;
        quality?: "S" | "A" | "B" | "C" | "D";
      }>;
    }>
  >([]);
  const [batchStatus, setBatchStatus] = useState<{
    total: number;
    done: number;
    cancelled?: boolean;
  } | null>(null);
  const batchCancelRef = useRef<boolean>(false);
  const blastCacheRef = useRef<Map<string, BlastResponse>>(new Map());
  const { setActiveTab, setPresetReversePair } = useWorkbench();
  const [showManualInputs, setShowManualInputs] = useState<boolean>(false);
  const [autoBlastExons, setAutoBlastExons] = useState<boolean>(true);
  const [autoBlastTopN, setAutoBlastTopN] = useState<number>(2);
  const normalizedSeq = useMemo(
    () =>
      sequence
        .replace(/^>.*/gm, "") // FASTA ヘッダを除去
        .replace(/\s+/g, "")
        .toUpperCase(),
    [sequence],
  );
  const storageKey = useMemo(
    () => (mode === "sequence" ? "seq_workbench_exon_panel_sequence" : "seq_workbench_exon_panel_exon"),
    [mode],
  );

  // 設定と配列をローカルストレージに自動保存/復元
  useEffect(() => {
    if (typeof window === "undefined") return;
    // F5 (reload) のときはローカルストレージをクリアして復元しない
    const isReload = (() => {
      try {
        const nav = performance.getEntriesByType("navigation")[0] as PerformanceNavigationTiming | undefined;
        if (nav && nav.type) return nav.type === "reload";
        const legacy = (performance as any).navigation;
        return legacy && legacy.type === 1;
      } catch {
        return false;
      }
    })();
    if (isReload) {
      try {
        window.localStorage.removeItem(storageKey);
      } catch {
        // ignore
      }
      return;
    }
    try {
      const raw = window.localStorage.getItem(storageKey);
      if (!raw) return;
      const saved = JSON.parse(raw) as Partial<{
        sequence: string;
        exonRanges: string;
        cdsRange: string;
        geneId: string;
        species: string;
        cdsSequenceInput: string;
        exonSequencesInput: string;
        fastaInput: string;
        sequencingMode: boolean;
        seqProductMin: number;
        seqProductMax: number;
        seqOverlap: number;
        seqMargin: number;
        customLocalDb: string;
        blastUseLocal: boolean;
        blastUseNcbi: boolean;
        ncbiOrganism: boolean;
        ncbiArabidopsis: boolean;
        ncbiCustomQuery: string;
      }>;
      if (saved.sequence) setSequence(saved.sequence);
      if (saved.exonRanges) setExonRanges(saved.exonRanges);
      if (saved.cdsRange) setCdsRange(saved.cdsRange);
      if (saved.geneId) setGeneId(saved.geneId);
      if (saved.species) setSpecies(saved.species);
      if (saved.cdsSequenceInput) setCdsSequenceInput(saved.cdsSequenceInput);
      if (saved.exonSequencesInput) setExonSequencesInput(saved.exonSequencesInput);
      if (saved.fastaInput) setFastaInput(saved.fastaInput);
      if (typeof saved.sequencingMode === "boolean") setSequencingMode(saved.sequencingMode);
      if (typeof saved.seqProductMin === "number") setSeqProductMin(saved.seqProductMin);
      if (typeof saved.seqProductMax === "number") setSeqProductMax(saved.seqProductMax);
      if (typeof saved.seqOverlap === "number") setSeqOverlap(saved.seqOverlap);
      if (typeof saved.seqMargin === "number") setSeqMargin(saved.seqMargin);
      if (typeof saved.customLocalDb === "string") setCustomLocalDb(saved.customLocalDb);
      if (typeof saved.blastUseLocal === "boolean") setBlastUseLocal(saved.blastUseLocal);
      if (typeof saved.blastUseNcbi === "boolean") setBlastUseNcbi(saved.blastUseNcbi);
      if (typeof saved.ncbiOrganism === "boolean") setNcbiOrganism(saved.ncbiOrganism);
      if (typeof saved.ncbiArabidopsis === "boolean") setNcbiArabidopsis(saved.ncbiArabidopsis);
      if (typeof saved.ncbiCustomQuery === "string") setNcbiCustomQuery(saved.ncbiCustomQuery);
    } catch {
      // ignore
    }
  }, [storageKey]);

  useEffect(() => {
    if (typeof window === "undefined") return;
    const payload = {
      sequence,
      exonRanges,
      cdsRange,
      geneId,
      species,
      cdsSequenceInput,
      exonSequencesInput,
      fastaInput,
      sequencingMode,
      seqProductMin,
      seqProductMax,
      seqOverlap,
      seqMargin,
      customLocalDb,
      blastUseLocal,
      blastUseNcbi,
      ncbiOrganism,
      ncbiArabidopsis,
      ncbiCustomQuery,
    };
    try {
      window.localStorage.setItem(storageKey, JSON.stringify(payload));
    } catch {
      // ignore
    }
  }, [
    storageKey,
    sequence,
    exonRanges,
    cdsRange,
    geneId,
    species,
    cdsSequenceInput,
    exonSequencesInput,
    fastaInput,
    sequencingMode,
    seqProductMin,
    seqProductMax,
    seqOverlap,
    seqMargin,
    customLocalDb,
    blastUseLocal,
    blastUseNcbi,
    ncbiOrganism,
    ncbiArabidopsis,
  ]);

  // --- Gene Loading Logic ---

  const handleLoadGene = async () => {
    // 1. Resolve DB
    let dbVal = normalizeLocalDbValue(loadGeneDb === "custom" ? loadGeneCustomDb : loadGeneDb);
    const gene = loadGeneName.trim();
    if (!gene) {
      showToast("遺伝子名を入力してください", "error");
      return;
    }
    setLoadGeneLoading(true);
    try {
      let hit;

      if (loadGeneDb === "auto") {
        // Auto Detect: Search all known local DBs
        const dbsToSearch = localDbOptions.map(o => o.value).filter(v => v !== "custom");
        if (dbsToSearch.length === 0) {
          throw new Error("検索可能なローカルDBが見つかりません");
        }

        // Batch queries
        const queries = dbsToSearch.map(db => ({ db, ids: [gene] }));
        const locs = await bioapiClient.fetchGeneLocations({ queries });

        hit = locs.find(x => x.found);
        if (!hit) {
          throw new Error(`遺伝子 ${gene} がどのDBにも見つかりませんでした`);
        }
        dbVal = hit.db; // Resolved DB
        showToast(`Found ${gene} in ${dbVal}`, "success");
      } else {
        if (!dbVal) {
          showToast("DBを選択してください", "error");
          setLoadGeneLoading(false);
          return;
        }
        const locs = await bioapiClient.fetchGeneLocations({
          queries: [{ db: dbVal, ids: [gene] }],
        });
        hit = locs.find((x) => x.found);
        if (!hit || !hit.chrom || hit.start == null || hit.end == null) {
          throw new Error(`遺伝子 ${gene} が見つかりませんでした (DB: ${dbVal})`);
        }
      }

      const { chrom, start, end } = hit;
      if (!chrom || start == null || end == null) throw new Error("座標情報が不足しています");

      // 2. Margin Calculation
      let margin = 2000;
      if (typeof loadGeneMargin === "number") {
        margin = loadGeneMargin;
      } else {
        // Auto Margin: max(2000, geneLen * 0.5), cap 10000
        const geneLen = (end - start) + 1;
        margin = Math.max(2000, Math.min(10000, Math.floor(geneLen * 0.5)));
      }

      const fetchStart = Math.max(1, start - margin);
      const fetchEnd = end + margin;

      // 2. Fetch Sequence
      // strand は指定せず、まずはプラス鎖として取得し、後で遺伝子の向きに合わせて反転は...しない。
      // FeatureSequenceView はあくまで「ゲノムのプラス鎖」として表示するのが通例だが、
      // ユーザーの利便性のため「その遺伝子の向き」で取得するオプションがあると良い。
      // ここでは仕様通り「遺伝子領域全体」を取得するが、strand は考慮せずプラス鎖として取得し、
      // gene model の strand 情報を使って FeatureSequenceView 側で表現...は難しい。
      // シンプルに「プラス鎖」で取得し、Exon 座標もそのままプラス鎖座標として扱う。
      // もし遺伝子がマイナス鎖なら、FeatureSequenceView 上で「Reverse」と表示されるべきだが、
      // 現在の FeatureSequenceView は配列そのものを表示する。
      // なので、常にプラス鎖で取得する。

      const seqRes = await bioapiClient.fetchLocalDbSequence({
        db: dbVal,
        entry: chrom,
        start: fetchStart,
        end: fetchEnd,
        strand: "plus",
      });

      const fullSeq = seqRes.sequence;
      const header = `${gene} on ${chrom}:${fetchStart}-${fetchEnd} (margin=${margin})`;

      // 3. Fetch Gene Model
      const modelRes = await bioapiClient.fetchRegionGeneModel({
        db: dbVal,
        entry: chrom,
        start: fetchStart,
        end: fetchEnd,
        gene_hint: gene,
      });

      // ターゲット遺伝子を探す
      const targetGene = modelRes.genes.find((g) => g.gene_name === gene || g.gene_id === gene) || modelRes.genes[0];

      if (!targetGene) {
        // モデルが見つからなくても配列はセットする
        setSequence(`>${header}\n${fullSeq}`);
        setGenomeHeader(header);
        setExonRanges("");
        setCdsRange("");
        setGeneId(gene);
        showToast("配列を取得しました（遺伝子モデルは見つかりませんでした）", "info");
        return;
      }

      // 4. Map coordinates to relative
      // seqRes.start should be fetchStart
      const offset = seqRes.start || 1;

      const mapRange = (gStart: number, gEnd: number) => {
        // genome coord -> relative (0-based) ? No, input for exonRanges is 1-based relative to sequence start.
        // sequence start is offset.
        // relative start = (gStart - offset) + 1
        const rStart = gStart - offset + 1;
        const rEnd = gEnd - offset + 1;
        return `${Math.max(1, rStart)}-${Math.min(fullSeq.length, rEnd)}`;
      };

      const exons = targetGene.exons.map((e) => mapRange(e.start, e.end)).join(",");
      // CDS は複数ある場合、最小〜最大を range とする（現在の UI仕様が単一 range のため）
      // あるいは sequencePrimerPanel の場合は使い方が違うかもしれないが、当面はこれで。
      let cdsString = "";
      if (targetGene.cds.length > 0) {
        const sortedCds = [...targetGene.cds].sort((a, b) => a.start - b.start);
        const cdsStart = sortedCds[0].start;
        const cdsEnd = sortedCds[sortedCds.length - 1].end;
        cdsString = mapRange(cdsStart, cdsEnd);
      }

      setSequence(`>${header}\n${fullSeq}`);
      setGenomeHeader(header);
      setExonRanges(exons);
      setCdsRange(cdsString);
      setGeneId(gene); // 内部 ID として保持

      showToast(`ロード完了: ${gene} (${targetGene.strand > 0 ? "+" : "-"})`, "success");

      // 5. Construct FASTA for input
      // targetGene.exons / .cds are in genomic coordinates.
      // fullSeq starts at fetchStart (which corresponds to index 0).
      // We need to extract substrings.
      const getSubseq = (gStart: number, gEnd: number) => {
        const s = gStart - fetchStart;
        const e = gEnd - fetchStart;
        if (s < 0 || e >= fullSeq.length) return ""; // margin外など
        return fullSeq.slice(s, e + 1);
      };

      const isReverse = targetGene.strand < 0;

      // cDNA: concat all exons
      const genomicExons = [...targetGene.exons].sort((a, b) => a.start - b.start);
      const splicedGenomic = genomicExons.map(e => getSubseq(e.start, e.end)).join("");
      const cdnaSeq = isReverse ? revComp(splicedGenomic) : splicedGenomic;

      // CDS
      let cdsSeq = "";
      if (targetGene.cds.length > 0) {
        const genomicCds = [...targetGene.cds].sort((a, b) => a.start - b.start);
        const splicedGenomicCds = genomicCds.map(c => getSubseq(c.start, c.end)).join("");
        cdsSeq = isReverse ? revComp(splicedGenomicCds) : splicedGenomicCds;
      }

      // Individual Exons
      const exonFastaEntries: string[] = [];
      if (isReverse) {
        // Minus strand: Exon 1 = RevComp(HighestCoordExon)
        const cnt = genomicExons.length;
        for (let i = 0; i < cnt; i++) {
          const gExon = genomicExons[cnt - 1 - i];
          const seq = revComp(getSubseq(gExon.start, gExon.end));
          exonFastaEntries.push(`>${gene}_exon${i + 1}\n${seq}`);
        }
      } else {
        // Plus strand
        genomicExons.forEach((gExon, i) => {
          const seq = getSubseq(gExon.start, gExon.end);
          exonFastaEntries.push(`>${gene}_exon${i + 1}\n${seq}`);
        });
      }

      const fastaBlocks = [
        `>${gene} cDNA`,
        cdnaSeq,
        `>${gene} CDS`,
        cdsSeq,
        ...exonFastaEntries
      ];

      // Filter out empty lines just in case
      const finalFasta = fastaBlocks.filter(Boolean).join("\n");
      setFastaInput(finalFasta);

    } catch (e) {
      showToast(e instanceof Error ? e.message : "ロードに失敗しました", "error");
    } finally {
      setLoadGeneLoading(false);
    }
  };

  const buildNcbiTargets = useCallback((): NCBITarget[] => {
    const targets: NCBITarget[] = [];
    if (ncbiOrganism) targets.push({ label: "crop", entrez_query: "", database: "nt" });
    if (ncbiArabidopsis)
      targets.push({ label: "arabidopsis", entrez_query: "Arabidopsis thaliana[Organism]", database: "nt" });
    const custom = ncbiCustomQuery.trim();
    if (custom) targets.push({ label: "custom", entrez_query: custom, database: "nt" });
    return targets;
  }, [ncbiOrganism, ncbiArabidopsis, ncbiCustomQuery]);
    const sidebarPrimers = useMemo(() => {
      const entries: Array<{
        kind?: "exon" | "cds";
        exonIndex: number;
        exonEndIndex?: number;
      subIndex?: number;
      left: string;
      right: string;
      product?: number;
      note?: string;
      quality?: "S" | "A" | "B" | "C" | "D";
      leftStart?: number;
      leftLen?: number;
      rightStart?: number;
      rightLen?: number;
    }> = [];
    exonPrimerPlans.forEach((p) => {
      const list =
        p.filteredCandidates?.slice(0, 2) ??
        (p.result
          ? p.result.candidates.slice(0, 2).map((c) => ({ pair: c, note: undefined }))
          : []);
      list.forEach((entry) => {
        const cand = entry.pair;
        const quality = (entry as any).quality as
          | "S"
          | "A"
          | "B"
          | "C"
          | "D"
          | undefined;
        entries.push({
          kind: p.kind,
          exonIndex: p.exonIndex,
          exonEndIndex: p.exonEndIndex,
          subIndex: p.subIndex,
          left: cand.left_sequence,
          right: cand.right_sequence,
          product: cand.product_size ?? undefined,
          note: entry.note,
          quality,
          leftStart: cand.left_start,
          leftLen: cand.left_length,
          rightStart: cand.right_start ?? cand.left_start,
          rightLen: cand.right_length,
        });
      });
    });
    return entries;
  }, [exonPrimerPlans]);

  const patternPrimers = useMemo(() => {
    const makeList = (pickIndex: number) => {
      const arr: typeof sidebarPrimers = [];
      exonPrimerPlans.forEach((p) => {
        const list =
          p.filteredCandidates ??
          (p.result ? p.result.candidates.map((c) => ({ pair: c, note: undefined })) : []);
        const entry = list[pickIndex];
        if (!entry) return;
        const cand = entry.pair;
        const quality = (entry as any).quality as
          | "S"
          | "A"
          | "B"
          | "C"
          | "D"
          | undefined;
        arr.push({
          kind: p.kind,
          exonIndex: p.exonIndex,
          exonEndIndex: p.exonEndIndex,
          subIndex: p.subIndex,
          left: cand.left_sequence,
          right: cand.right_sequence,
          product: cand.product_size ?? undefined,
          note: entry.note,
          quality,
          leftStart: cand.left_start,
          leftLen: cand.left_length,
          rightStart: cand.right_start ?? cand.left_start,
          rightLen: cand.right_length,
        });
      });
      return arr;
    };
    return {
      pattern1: makeList(0),
      pattern2: makeList(1),
    };
  }, [exonPrimerPlans]);
  const currentPatternPrimers = patternPrimers[primerPatternTab] ?? [];
  const topCandidateSummary = useMemo(() => {
    return exonPrimerPlans
      .filter(
        (p) =>
          (p.filteredCandidates && p.filteredCandidates.length > 0) ||
          (p.result && p.result.candidates.length > 0),
      )
      .map((p) => {
        const cand =
          p.filteredCandidates?.[0]?.pair ??
          (p.result ? p.result.candidates[0] : null);
        if (!cand) return null;
        const firstFiltered = p.filteredCandidates?.[0] as
          | { quality?: "S" | "A" | "B" | "C" | "D" }
          | undefined;
        return {
          kind: p.kind ?? "exon",
          exonIndex: p.exonIndex,
          range: p.range,
          left: cand.left_sequence,
          right: cand.right_sequence,
          product: cand.product_size,
          penalty: cand.pair_penalty,
          quality: firstFiltered?.quality,
        };
      })
      .filter((v): v is NonNullable<typeof v> => !!v);
  }, [exonPrimerPlans]);

  const parsedExons = useMemo(() => parseRanges(exonRanges), [exonRanges]);
  const parsedCds = useMemo(() => parseRanges(cdsRange), [cdsRange]);

  const featureViewExons = useMemo<Range1Based[]>(
    () => parsedExons.map(([start, end]) => ({ start, end })),
    [parsedExons],
  );
  const featureViewCds = useMemo<Range1Based[]>(
    () => parsedCds.map(([start, end]) => ({ start, end })),
    [parsedCds],
  );

  const featureViewPrimerRanges = useMemo<PrimerRange1Based[]>(() => {
    const out: PrimerRange1Based[] = [];
    currentPatternPrimers.forEach((p) => {
      if (p.leftStart != null && p.leftLen) {
        out.push({ start: p.leftStart, end: p.leftStart + p.leftLen - 1, kind: "left" });
      }
      if (p.rightStart != null && p.rightLen) {
        out.push({ start: p.rightStart, end: p.rightStart + p.rightLen - 1, kind: "right" });
      }
    });
    return out;
  }, [currentPatternPrimers]);

  const coverageStats = useMemo(() => {
    if (!normalizedSeq) return null;
    const primers = currentPatternPrimers;
    if (!primers.length) return null;
    const n = normalizedSeq.length;
    const covered = new Array<boolean>(n).fill(false);

    const clampRange = (start: number, end: number): [number, number] => {
      let s = Math.max(1, start);
      let e = Math.min(n, end);
      if (e < s) [s, e] = [e, s];
      return [s, e];
    };

    const computeAmpRange = (p: {
      leftStart?: number;
      leftLen?: number;
      rightStart?: number;
      rightLen?: number;
      product?: number;
    }): [number, number] => {
      const leftStart = p.leftStart ?? 1;
      const leftEnd = p.leftLen ? leftStart + p.leftLen - 1 : leftStart;
      const rightStart = p.rightStart ?? leftStart;
      const rightEnd = p.rightLen ? rightStart + p.rightLen - 1 : rightStart;
      let regionStart = Math.min(leftStart, rightStart);
      let regionEnd = Math.max(leftEnd, rightEnd);
      if (p.product && p.product > 0) {
        regionEnd = regionStart + p.product - 1;
      }
      return clampRange(regionStart, regionEnd);
    };

    primers.forEach((p) => {
      const [s, e] = computeAmpRange(p);
      for (let i = s - 1; i < e; i += 1) {
        covered[i] = true;
      }
    });

    const countInRanges = (ranges: Array<[number, number]>) => {
      let total = 0;
      let hit = 0;
      ranges.forEach(([sRaw, eRaw]) => {
        const [s, e] = clampRange(sRaw, eRaw);
        for (let pos = s; pos <= e; pos += 1) {
          total += 1;
          if (covered[pos - 1]) hit += 1;
        }
      });
      return { total, hit };
    };

    const exon = parsedExons.length ? countInRanges(parsedExons) : null;
    const cds = parsedCds.length ? countInRanges(parsedCds) : null;
    const genomeTotal = n;
    const genomeHit = covered.reduce((acc, v) => (v ? acc + 1 : acc), 0);

    return {
      exon,
      cds,
      genome: { total: genomeTotal, hit: genomeHit },
    };
  }, [normalizedSeq, currentPatternPrimers, parsedExons, parsedCds]);

  const labelForDb = (path: string) =>
    labelForDbPath(path, localDbOptions);

  const toggleLocalDb = (path: string) => {
    setSelectedLocalDbs((prev) =>
      prev.includes(path) ? prev.filter((p) => p !== path) : [...prev, path],
    );
  };

  const runBlastForSequence = async (seq: string): Promise<BlastResponse> => {
    const key = JSON.stringify({
      seq,
      dbs: effectiveLocalDbs,
      useLocal: blastUseLocal,
      useNcbi: blastUseNcbi,
      local_mode: localMode,
      max_target_seqs: blastMaxHits,
      max_hsps: blastMaxHsps,
      task: blastTask,
      evalue: blastEvalue,
    });
    const cached = blastCacheRef.current.get(key);
    if (cached) {
      return cached;
    }

    const localHits: BlastHit[] = [];
    if (blastUseLocal) {
      const batch = await bioapiClient.runBlastBatchLocal({
        sequences: [seq],
        dbs: effectiveLocalDbs,
        local_mode: localMode,
        task: blastTask,
        evalue: blastEvalue,
        max_target_seqs: blastMaxHits,
        max_hsps: blastMaxHsps ?? undefined,
        num_threads: blastNumThreads ?? undefined,
      });
      const merged = batch.results?.[0] ?? null;
      localHits.push(...relabelLocalBlastHits(merged?.hits ?? [], effectiveLocalDbs, localDbOptions));
    }

    let ncbiHits: BlastHit[] = [];
    if (blastUseNcbi) {
      const targets = buildNcbiTargets();
      if (targets.length === 0) {
        throw new Error("NCBI を使う場合、ターゲット種またはクエリを 1 つ以上選択してください。");
      }
      const res = await bioapiClient.runBlastMulti({
        sequence: seq,
        db: effectiveLocalDbs[0] ?? "",
        max_target_seqs: blastMaxHits,
        max_hsps: blastMaxHsps ?? undefined,
        task: "blastn",
        evalue: blastEvalue,
        num_threads: blastNumThreads ?? undefined,
        backends: ["ncbi"],
        ncbi_targets: targets,
        ncbi_database: "nt",
      });
      ncbiHits = (res.hits ?? []).map((h) => ({ ...h, source: h.source || "ncbi" }));
    }

    const hits = [...localHits, ...ncbiHits];
    const res: BlastResponse = { num_hits: hits.length, hits };
    blastCacheRef.current.set(key, res);
    return res;
  };

  const revComp = (s: string) =>
    s
      .split("")
      .reverse()
      .map((c) => {
        switch (c) {
          case "A":
            return "T";
          case "T":
            return "A";
          case "G":
            return "C";
          case "C":
            return "G";
          default:
            return c;
        }
      })
      .join("");

  const fetchStructure = async () => {
    if (!geneId.trim()) {
      setStructureError("Gene ID を入力してください。");
      return;
    }
    setStructureLoading(true);
    setStructureError(null);
    try {
      const data = await bioapiClient.fetchGeneStructure(geneId.trim(), species.trim() || undefined);
      setSequence(data.sequence);
      if (data.exons?.length) {
        setExonRanges(data.exons.map((r) => `${r.start}-${r.end}`).join(","));
      }
      if (data.cds?.length) {
        setCdsRange(data.cds.map((r) => `${r.start}-${r.end}`).join(","));
        // 代表CDSを自動選択
        const cds0 = data.cds[0];
        setSelectionRange({ start: cds0.start, end: cds0.end });
        setSubStart(cds0.start);
        setSubEnd(cds0.end);
      }
    } catch (e) {
      setStructureError(e instanceof Error ? e.message : "遺伝子構造の取得に失敗しました。");
    } finally {
      setStructureLoading(false);
    }
  };

  const guessGeneIdByBlast = async () => {
    if (!normalizedSeq) {
      setStructureError("まず配列を入力してください。");
      return;
    }
    if (!blastUseLocal && !blastUseNcbi) {
      setStructureError("少なくとも1つは BLAST 実行先を選択してください。");
      return;
    }
    if (blastUseLocal && effectiveLocalDbs.length === 0) {
      setStructureError("ローカル BLAST+ を使う場合は DB を選択してください。");
      return;
    }
    setStructureLoading(true);
    setStructureError(null);
    setBlastJobId(null);
    setBlastJobInfo(null);
    try {
      const hits: BlastHit[] = [];
      if (blastUseLocal) {
        const { result: batch } = await runBlastBatchLocalJob(
          {
            sequences: [normalizedSeq],
            dbs: effectiveLocalDbs,
            local_mode: localMode,
            task: blastTask,
            evalue: blastEvalue,
            max_target_seqs: blastMaxHits,
            max_hsps: blastMaxHsps ?? undefined,
            num_threads: blastNumThreads ?? undefined,
            engine: "blast",
          },
          {
            onCreated: (id) => setBlastJobId(id),
            onUpdate: (info) => setBlastJobInfo(info),
            intervalMs: 900,
          },
        );
        const merged = batch.results?.[0] ?? null;
        hits.push(...relabelLocalBlastHits(merged?.hits ?? [], effectiveLocalDbs, localDbOptions));
      }

      if (blastUseNcbi) {
        const targets = buildNcbiTargets();
        if (targets.length === 0) {
          throw new Error("NCBI を使う場合、ターゲット種またはクエリを 1 つ以上選択してください。");
        }
        const resNcbi = await bioapiClient.runBlastMulti({
          sequence: normalizedSeq,
          db: effectiveLocalDbs[0] ?? "",
          max_target_seqs: blastMaxHits,
          max_hsps: blastMaxHsps ?? undefined,
          task: "blastn",
          evalue: blastEvalue,
          num_threads: blastNumThreads ?? undefined,
          backends: ["ncbi"],
          ncbi_targets: targets,
          ncbi_database: "nt",
        });
        hits.push(...((resNcbi.hits ?? []).map((h) => ({ ...h, source: h.source || "ncbi" })) as BlastHit[]));
      }

      hits.sort((a, b) => {
        if (a.bitscore !== b.bitscore) return b.bitscore - a.bitscore;
        return b.pident - a.pident;
      });
      if (hits.length === 0) {
        setStructureError("BLASTでヒットが見つかりませんでした。");
        return;
      }
      const top = hits[0].sseqid.split(/\s+/)[0];
      setGeneGuess(top);
      setGeneId(top);
    } catch (e) {
      setStructureError(e instanceof Error ? e.message : "BLAST で Gene 候補取得に失敗しました。");
    } finally {
      setStructureLoading(false);
      setBlastJobId(null);
      setBlastJobInfo(null);
    }
  };

  const FEATURE_VIEW_BLOCK_LEN = 150;

  const scrollToRegion = (start: number, end: number) => {
    const container = fastaHighlightRef.current;
    if (!container || !normalizedSeq) return;
    const centerPos = Math.max(
      1,
      Math.min(normalizedSeq.length, Math.floor((start + end) / 2)),
    );
    const blockStart = Math.floor((centerPos - 1) / FEATURE_VIEW_BLOCK_LEN) * FEATURE_VIEW_BLOCK_LEN + 1;
    const target = container.querySelector<HTMLElement>(`#feature-seq-${blockStart}`);
    if (!target) return;
    target.scrollIntoView({ block: "center", behavior: "smooth" });
    const ts = Date.now();
    setGuideRegion({ start, end, ts });
    window.setTimeout(() => {
      setGuideRegion((prev) => {
        if (!prev || prev.ts !== ts) return prev;
        return null;
      });
    }, 1500);
  };

  const handleDesign = async () => {
    if (!normalizedSeq) {
      setDesignError("まず配列を入力してください。");
      return;
    }
    const tStart = selectionRange
      ? selectionRange.start
      : subStart && subEnd
        ? Math.max(1, Math.min(subStart, subEnd))
        : null;
    const tLen = selectionRange
      ? Math.max(1, selectionRange.end - selectionRange.start + 1)
      : subStart && subEnd
        ? Math.max(1, Math.abs(subEnd - subStart) + 1)
        : null;

    setDesignLoading(true);
    setDesignError(null);
    setDesignResult(null);
    setSelectedPair(null);
    try {
      let res = await bioapiClient.designPrimers({
        sequence: normalizedSeq,
        num_return: numReturn,
        product_size_range: productSizeRange || null,
        target_start: tStart,
        target_length: tLen,
        opt_tm: pcrOptTm,
        min_tm: pcrMinTm,
        max_tm: pcrMaxTm,
        primer_min_size: pcrMinSize,
        primer_opt_size: pcrOptSize,
        primer_max_size: pcrMaxSize,
        primer_min_gc: pcrMinGc,
        primer_max_gc: pcrMaxGc,
        primer_salt_monovalent: pcrSaltMonovalent,
        primer_dna_conc: pcrDnaConc,
      });
      // 1回目で候補が 0 の場合は、product size 制約を緩めて再トライ（CDS/エキソン増幅タブのみ）
      if (mode === "exon" && !res.candidates.length) {
        const relaxed = await bioapiClient.designPrimers({
          sequence: normalizedSeq,
          num_return: numReturn,
          product_size_range: null,
          target_start: tStart,
          target_length: tLen,
          opt_tm: pcrOptTm,
          min_tm: pcrMinTm,
          max_tm: pcrMaxTm,
          primer_min_size: pcrMinSize,
          primer_opt_size: pcrOptSize,
          primer_max_size: pcrMaxSize,
          primer_min_gc: pcrMinGc,
          primer_max_gc: pcrMaxGc,
          primer_salt_monovalent: pcrSaltMonovalent,
          primer_dna_conc: pcrDnaConc,
        });
        if (relaxed.candidates.length) {
          res = relaxed;
        }
      }
      setDesignResult(res);
      setSelectedPair(res.candidates.length > 0 ? 0 : null);
    } catch (e) {
      setDesignError(e instanceof Error ? e.message : "プライマー設計に失敗しました。");
    } finally {
      setDesignLoading(false);
    }
  };

  const renderAmplicon = () => {
    if (!designResult || selectedPair === null) return null;
    const pair = designResult.candidates[selectedPair];
    const ampStart = Math.max(1, Math.min(pair.left_start, pair.right_start ?? pair.left_start));
    const ampEnd = pair.product_size ? ampStart + pair.product_size - 1 : Math.max(
      pair.left_start + pair.left_length - 1,
      (pair.right_start ?? pair.left_start) + pair.right_length - 1,
    );
    const prefix = normalizedSeq.slice(0, ampStart - 1);
    const amp = normalizedSeq.slice(ampStart - 1, ampEnd);
    const suffix = normalizedSeq.slice(ampEnd);
    return (
      <div className="amplicon-viewer">
        <p className="amplicon-title">
          選択ペアの増幅領域: {ampStart}–{ampEnd} (長さ {ampEnd - ampStart + 1} bp)
        </p>
        <div className="amplicon-seq">
          <span>{prefix}</span>
          <span className="amplicon-highlight">{amp}</span>
          <span>{suffix}</span>
        </div>
      </div>
    );
  };

  const cancelBlastJob = async () => {
    if (!blastJobId) return;
    try {
      const updated = await bioapiClient.cancelJob(blastJobId);
      setBlastJobInfo(updated);
    } catch (e) {
      setBlastError(e instanceof Error ? e.message : "ジョブのキャンセルに失敗しました。");
    }
  };

  const runPrimerBlast = async () => {
    if (!designResult || selectedPairIndices.length === 0) {
      setBlastError("BLASTするペアを選択してください。");
      return;
    }
    if (!blastUseLocal && !blastUseNcbi) {
      setBlastError("少なくとも1つはBLAST実行先を選んでください。");
      return;
    }
    if (blastUseLocal && effectiveLocalDbs.length === 0) {
      setBlastError("ローカルDBを使う場合は DB を選択してください。");
      return;
    }
    setBlastLoading(true);
    setBlastError(null);
    setBlastLeft(null);
    setBlastRight(null);
    setPrimerBlastResults([]);
    setBlastJobId(null);
    setBlastJobInfo(null);
    try {
      const selectedPairs = selectedPairIndices
        .map((idx) => ({ idx, pair: designResult.candidates[idx] }))
        .filter((x) => x.pair?.left_sequence && x.pair?.right_sequence);

      const uniqPrimerSeqs: string[] = [];
      const idxBySeq = new Map<string, number>();
      const indexOf = (raw: string): string => {
        const key = raw.replace(/\s+/g, "").toUpperCase();
        if (!key) return "";
        if (idxBySeq.has(key)) return key;
        idxBySeq.set(key, uniqPrimerSeqs.length);
        uniqPrimerSeqs.push(key);
        return key;
      };

      const pairRefs = selectedPairs.map(({ idx, pair }) => ({
        pairIndex: idx,
        left: indexOf(pair.left_sequence),
        right: indexOf(pair.right_sequence),
      }));

      const localMap = new Map<string, BlastResponse>();
      if (blastUseLocal) {
        const { result: batch } = await runBlastBatchLocalJob(
          {
            sequences: uniqPrimerSeqs,
            dbs: effectiveLocalDbs,
            local_mode: localMode,
            task: blastTask,
            evalue: blastEvalue,
            max_target_seqs: blastMaxHits,
            max_hsps: blastMaxHsps ?? undefined,
            num_threads: blastNumThreads ?? undefined,
            engine: "blast",
          },
          {
            onCreated: (id) => setBlastJobId(id),
            onUpdate: (info) => setBlastJobInfo(info),
            intervalMs: 900,
          },
        );
        if (!batch.results || batch.results.length !== uniqPrimerSeqs.length) {
          throw new Error("ローカル BLAST の結果件数が期待値と一致しませんでした。");
        }
        for (let i = 0; i < uniqPrimerSeqs.length; i += 1) {
          const seq = uniqPrimerSeqs[i];
          const merged = batch.results[i];
          const hits = relabelLocalBlastHits(merged?.hits ?? [], effectiveLocalDbs, localDbOptions);
          localMap.set(seq, { num_hits: hits.length, hits });
        }
      }

      const ncbiMap = new Map<string, BlastResponse>();
      if (blastUseNcbi) {
        const targets = buildNcbiTargets();
        if (targets.length === 0) {
          throw new Error("NCBI を使う場合、ターゲット種またはクエリを 1 つ以上選択してください。");
        }
        for (const seq of uniqPrimerSeqs) {
          const res = await bioapiClient.runBlastMulti({
            sequence: seq,
            db: effectiveLocalDbs[0] ?? "",
            max_target_seqs: blastMaxHits,
            max_hsps: blastMaxHsps ?? undefined,
            task: "blastn",
            evalue: blastEvalue,
            num_threads: blastNumThreads ?? undefined,
            backends: ["ncbi"],
            ncbi_targets: targets,
            ncbi_database: "nt",
          });
          const hits = (res.hits ?? []).map((h) => ({ ...h, source: h.source || "ncbi" }));
          ncbiMap.set(seq, { num_hits: hits.length, hits });
        }
      }

      const results = pairRefs.map((ref) => {
        const leftLocal = localMap.get(ref.left) ?? { num_hits: 0, hits: [] };
        const rightLocal = localMap.get(ref.right) ?? { num_hits: 0, hits: [] };
        const leftNcbi = ncbiMap.get(ref.left) ?? { num_hits: 0, hits: [] };
        const rightNcbi = ncbiMap.get(ref.right) ?? { num_hits: 0, hits: [] };
        const leftHits = [...(leftLocal.hits ?? []), ...(leftNcbi.hits ?? [])];
        const rightHits = [...(rightLocal.hits ?? []), ...(rightNcbi.hits ?? [])];
        return {
          pairIndex: ref.pairIndex,
          left: { num_hits: leftHits.length, hits: leftHits },
          right: { num_hits: rightHits.length, hits: rightHits },
        };
      });

      setPrimerBlastResults(results);
    } catch (e) {
      setBlastError(e instanceof Error ? e.message : "Primer-BLAST 実行中にエラーが発生しました");
    } finally {
      setBlastLoading(false);
    }
  };

  const designPrimersForExons = async () => {
    if (!normalizedSeq) {
      setDesignError("配列が読み込まれていません。");
      return;
    }
    if (parsedExons.length === 0) {
      setDesignError("エキソン情報がありません。FASTA ヘッダーまたは入力欄を確認してください。");
      return;
    }
    setDesignError(null);

    const useSequencingLogic = mode === "sequence" ? true : sequencingMode;

    // Sequencing Mode Logic
    let plansToDesign: Array<{
      kind?: "exon" | "cds";
      exonIndex: number;
      exonEndIndex?: number;
      subIndex?: number;
      range: [number, number];
      result?: PrimerDesignResponse;
      error?: string;
      loading?: boolean;
      filteredCandidates?: Array<{
        pair: PrimerPair;
        left?: BlastResponse;
        right?: BlastResponse;
        note?: string;
        ampliconCount?: number;
        quality?: "S" | "A" | "B" | "C" | "D";
      }>;
    }> = [];

    const addTiledRanges = (start: number, end: number, exonIndex: number) => {
      let currentStart = start;
      let subIdx = 1;
      while (currentStart < end) {
        let currentEnd = currentStart + seqProductMax - 1;
        if (currentEnd > end) currentEnd = end;

        if (currentEnd - currentStart + 1 < seqProductMin) {
          currentStart = Math.max(start, currentEnd - seqProductMin + 1);
        }

        plansToDesign.push({
          kind: "exon",
          exonIndex,
          subIndex: subIdx,
          range: [currentStart, currentEnd],
          loading: true,
        });

        if (currentEnd >= end) break;
        currentStart = currentEnd - seqOverlap + 1;
        subIdx += 1;
      }
    };

    if (useSequencingLogic) {
      if (mode === "sequence") {
        // シーケンスプライマー用:
        // - まずは 600–800bp の産物で、複数エキソンを 1 ペアにまとめられるところはまとめる。
        // - それでも 1 エキソン (+margin) が長すぎる場合は、そのエキソンだけをタイル分割。
        const n = parsedExons.length;
        let i = 0;

        while (i < n) {
          const [sRaw, eRaw] = parsedExons[i];
          const startWithMargin = Math.max(1, sRaw - seqMargin);
          const endWithMargin = Math.min(normalizedSeq.length, eRaw + seqMargin);
          const spanLen = endWithMargin - startWithMargin + 1;

          if (spanLen > seqProductMax) {
            // このエキソン単体が長すぎるので、エキソンごとにタイル分割
            addTiledRanges(startWithMargin, endWithMargin, i + 1);
            i += 1;
            continue;
          }

          let blockStart = startWithMargin;
          let blockEnd = endWithMargin;
          const blockExonIndex = i + 1;
          let blockExonEndIndex = i + 1;
          let j = i + 1;

          while (j < n) {
            const [nextS, nextE] = parsedExons[j];
            const nextStart = Math.max(1, nextS - seqMargin);
            const nextEnd = Math.min(normalizedSeq.length, nextE + seqMargin);
            const newSpan = nextEnd - blockStart + 1;
            if (newSpan <= seqProductMax) {
              blockEnd = nextEnd;
              blockExonEndIndex = j + 1;
              j += 1;
            } else {
              break;
            }
          }

          const blockLen = blockEnd - blockStart + 1;
          if (blockLen <= seqProductMax) {
            plansToDesign.push({
              kind: "exon",
              exonIndex: blockExonIndex,
              exonEndIndex: blockExonEndIndex,
              range: [blockStart, blockEnd],
              loading: true,
            });
          } else {
            addTiledRanges(blockStart, blockEnd, blockExonIndex);
          }

          i = j;
        }
      } else {
        // CDS/エキソン増幅タブで Sequencing Mode が ON の場合:
        // 600〜800bp の産物で、できるだけ多くのエキソンを 1 ペアでカバーするようにまとめる。
        const n = parsedExons.length;
        let i = 0;

        while (i < n) {
          const [sRaw, eRaw] = parsedExons[i];
          const startWithMargin = Math.max(1, sRaw - seqMargin);
          const endWithMargin = Math.min(normalizedSeq.length, eRaw + seqMargin);
          const spanLen = endWithMargin - startWithMargin + 1;

          // 単一エキソン（+マージン）が max を超える場合は、そのエキソンだけをタイル分割
          if (spanLen > seqProductMax) {
            addTiledRanges(startWithMargin, endWithMargin, i + 1);
            i += 1;
            continue;
          }

          // まず現在のエキソンを含むブロックを開始
          let blockStart = startWithMargin;
          let blockEnd = endWithMargin;
          const blockExonIndex = i + 1;
          let blockExonEndIndex = i + 1;
          let j = i + 1;

          // 次のエキソンをできるだけ多く同じブロックに詰め込む（産物長 <= max）
          while (j < n) {
            const [nextS, nextE] = parsedExons[j];
            const nextStart = Math.max(1, nextS - seqMargin);
            const nextEnd = Math.min(normalizedSeq.length, nextE + seqMargin);
            const newSpan = nextEnd - blockStart + 1;
            if (newSpan <= seqProductMax) {
              blockEnd = nextEnd;
              blockExonEndIndex = j + 1;
              j += 1;
            } else {
              break;
            }
          }

          const blockLen = blockEnd - blockStart + 1;
          if (blockLen <= seqProductMax) {
            const ampStart = blockStart;
            const ampEnd = blockEnd;
            plansToDesign.push({
              kind: "exon",
              exonIndex: blockExonIndex,
              exonEndIndex: blockExonEndIndex,
              range: [ampStart, ampEnd],
              loading: true,
            });
          } else {
            // 理論上ここには来ないはずだが、安全のためタイル分割
            addTiledRanges(blockStart, blockEnd, blockExonIndex);
          }

          i = j;
        }
      }
    } else {
      plansToDesign = parsedExons.map((r, idx) => ({
        kind: "exon",
        exonIndex: idx + 1,
        range: r,
        loading: true,
      }));
    }

    setExonPrimerPlans(plansToDesign);
    batchCancelRef.current = false;
    setBatchStatus(
      plansToDesign.length
        ? { total: plansToDesign.length, done: 0, cancelled: false }
        : null,
    );

    const hw = typeof navigator !== "undefined" ? navigator.hardwareConcurrency || 12 : 12;
    const baseParallel = Math.max(1, Math.min(12, Math.floor(hw / 2) || 1));
    const parallel = Math.max(
      1,
      Math.min(
        plansToDesign.length,
        autoBlastExons ? Math.min(4, baseParallel) : baseParallel,
      ),
    );

    const runPlan = async (
      plan: (typeof plansToDesign)[number],
    ): Promise<(typeof plansToDesign)[number]> => {
      let [s, e] = plan.range;
      try {
        const len = Math.max(1, e - s + 1);
        let productRange = "";

        if (useSequencingLogic) {
          productRange = `${seqProductMin}-${seqProductMax}`;
        } else {
          const lower = Math.max(80, len);
          const upper = Math.min(normalizedSeq.length, lower + 150);
          productRange = `${lower}-${upper}`;
        }

        let res = await bioapiClient.designPrimers({
          sequence: normalizedSeq,
          num_return: Math.min(numReturn, 5),
          product_size_range: productRange,
          target_start: s,
          target_length: len,
          opt_tm: pcrOptTm,
          min_tm: pcrMinTm,
          max_tm: pcrMaxTm,
          primer_min_size: pcrMinSize,
          primer_opt_size: pcrOptSize,
          primer_max_size: pcrMaxSize,
          primer_min_gc: pcrMinGc,
          primer_max_gc: pcrMaxGc,
          primer_salt_monovalent: pcrSaltMonovalent,
          primer_dna_conc: pcrDnaConc,
        });
        // 1回目で候補が 0 の場合:
        // - Sequencing ロジック有効時（シーケンスプライマー / Sequencing Mode）:
        //   エキソン全体（+margin）内で target_start / target_length を動かしながら再トライ
        //   ※ PRIMER_PRODUCT_SIZE_RANGE（例: 600–800bp）は維持する
        // - 通常の CDS/エキソン増幅タブ:
        //   product size 制約を緩めて再トライ
        if (useSequencingLogic && !res.candidates.length) {
          const rangesForKind = plan.kind === "cds" ? parsedCds : parsedExons;
          const span = rangesForKind[plan.exonIndex - 1];
          if (span) {
            const [spanStart, spanEnd] = span;
            const allowStart = Math.max(1, spanStart - seqMargin);
            const allowEnd = Math.min(normalizedSeq.length, spanEnd + seqMargin);
            const maxWindow = seqProductMax;
            const minWindow = seqProductMin;

            // target_length の候補（元の長さから、600bp 付近まで段階的に短くする）
            const lengthsToTry: number[] = [];
            const initialLen = Math.min(len, maxWindow);
            const minTargetLen = Math.max(
              Math.min(len, minWindow),
              Math.floor(minWindow * 0.6),
            );
            lengthsToTry.push(initialLen);
            for (let t = initialLen - 50; t >= minTargetLen; t -= 50) {
              if (!lengthsToTry.includes(t)) lengthsToTry.push(t);
            }

            let best: PrimerDesignResponse | null = null;
            let bestStart = s;
            let bestLen = len;
            let trialCount = 0;
            const spanLength = allowEnd - allowStart + 1;
            const step = 30;
            const approxPerLen = Math.max(1, Math.ceil(spanLength / step) * 2);
            const maxTrials = Math.min(800, approxPerLen * Math.max(1, lengthsToTry.length));

            outer: for (const tLen of lengthsToTry) {
              const maxStart = allowEnd - tLen + 1;
              const minStart = allowStart;
              if (maxStart < minStart) {
                continue;
              }
              const centerStart = Math.min(Math.max(s, minStart), maxStart);

              for (let offset = 0; ; offset += step) {
                const candidateStarts =
                  offset === 0 ? [centerStart] : [centerStart - offset, centerStart + offset];
                for (const newStart of candidateStarts) {
                  if (newStart < minStart || newStart > maxStart) continue;
                  const newEnd = newStart + tLen - 1;
                  if (newEnd > allowEnd) continue;

                  trialCount += 1;
                  if (trialCount > maxTrials) {
                    break outer;
                  }

                  const trial = await bioapiClient.designPrimers({
                    sequence: normalizedSeq,
                    num_return: Math.min(numReturn, 5),
                    product_size_range: productRange,
                    target_start: newStart,
                    target_length: tLen,
                    opt_tm: pcrOptTm,
                    min_tm: pcrMinTm,
                    max_tm: pcrMaxTm,
                    primer_min_size: pcrMinSize,
                    primer_opt_size: pcrOptSize,
                    primer_max_size: pcrMaxSize,
                    primer_min_gc: pcrMinGc,
                    primer_max_gc: pcrMaxGc,
                    primer_salt_monovalent: pcrSaltMonovalent,
                    primer_dna_conc: pcrDnaConc,
                  });
                  if (trial.candidates.length) {
                    best = trial;
                    bestStart = newStart;
                    bestLen = tLen;
                    break outer;
                  }
                }
                if (centerStart - offset <= minStart && centerStart + offset >= maxStart) {
                  // 許容範囲を一通り走査し終えた
                  break;
                }
              }
            }

            if (best) {
              res = best;
              s = bestStart;
              e = bestStart + bestLen - 1;
            }
          }
        } else if (!useSequencingLogic && mode === "exon" && !res.candidates.length) {
          const relaxed = await bioapiClient.designPrimers({
            sequence: normalizedSeq,
            num_return: Math.min(numReturn, 5),
            product_size_range: null,
            target_start: s,
            target_length: len,
            opt_tm: pcrOptTm,
            min_tm: pcrMinTm,
            max_tm: pcrMaxTm,
            primer_min_size: pcrMinSize,
            primer_opt_size: pcrOptSize,
            primer_max_size: pcrMaxSize,
            primer_min_gc: pcrMinGc,
            primer_max_gc: pcrMaxGc,
            primer_salt_monovalent: pcrSaltMonovalent,
            primer_dna_conc: pcrDnaConc,
          });
          if (relaxed.candidates.length) {
            res = relaxed;
          }
        }
        let filteredCandidates:
          | Array<{
            pair: PrimerPair;
            left: BlastResponse;
            right: BlastResponse;
            note?: string;
            ampliconCount?: number;
            quality?: "S" | "A" | "B" | "C" | "D";
          }>
          | undefined;

        if (autoBlastExons) {
          const scored: Array<{
            pair: PrimerPair;
            left: BlastResponse;
            right: BlastResponse;
            totalHits: number;
            amplicons: number;
          }> = [];
          const candidatesToEval = res.candidates.slice(
            0,
            Math.max(autoBlastTopN * 2, autoBlastTopN, 3),
          );

          const gradeFromAmp = (n: number): "S" | "A" | "B" | "C" | "D" => {
            if (n === 1) return "S";
            if (n === 0) return "D";
            if (n === 2) return "C";
            if (n === 3) return "B";
            return "D";
          };

          if (blastUseLocal && effectiveLocalDbs.length > 0) {
            const idxBySeq = new Map<string, number>();
            const uniq: string[] = [];
            const normalize = (s: string) => s.replace(/\s+/g, "").toUpperCase();
            const indexOf = (s: string) => {
              const key = normalize(s);
              const existing = idxBySeq.get(key);
              if (existing != null) return existing;
              const idx = uniq.length;
              idxBySeq.set(key, idx);
              uniq.push(key);
              return idx;
            };

            candidatesToEval.forEach((p) => {
              indexOf(p.left_sequence);
              indexOf(p.right_sequence);
            });

            if (uniq.length > 0) {
              const batch = await bioapiClient.runBlastBatchLocal({
                sequences: uniq,
                dbs: effectiveLocalDbs,
                local_mode: localMode,
                task: blastTask,
                evalue: blastEvalue,
                max_target_seqs: blastMaxHits,
                max_hsps: blastMaxHsps ?? undefined,
                num_threads: blastNumThreads ?? undefined,
              });
              if (!batch.results || batch.results.length !== uniq.length) {
                throw new Error("run_batch_local の結果件数が期待値と一致しませんでした。");
              }

              const perSeq = batch.results.map((r) => ({
                ...r,
                hits: relabelLocalBlastHits(r.hits ?? [], effectiveLocalDbs, localDbOptions),
              }));

              for (const pair of candidatesToEval) {
                const left = perSeq[indexOf(pair.left_sequence)];
                const right = perSeq[indexOf(pair.right_sequence)];
                const totalHits = countLocalHits(left) + countLocalHits(right);
                const { amplicons } = computePrimerAmplicons(left, right);
                scored.push({ pair, left, right, totalHits, amplicons: amplicons.length });
              }
            }
          } else {
            for (const pair of candidatesToEval) {
              const [left, right] = await Promise.all([
                runBlastForSequence(pair.left_sequence),
                runBlastForSequence(pair.right_sequence),
              ]);
              const totalHits = countLocalHits(left) + countLocalHits(right);
              const { amplicons } = computePrimerAmplicons(left, right);
              scored.push({ pair, left, right, totalHits, amplicons: amplicons.length });
            }
          }

          // まずは「予測PCR産物が 1 本」のペアを優先
          const strict = scored.filter((s) => s.amplicons === 1);
          let picked = strict;
          let note: string | undefined = "ローカルで予測PCR産物が1本のみ";

          // 1本が無ければ、1〜2本のペアを暫定採用
          if (picked.length === 0) {
            picked = scored.filter((s) => s.amplicons >= 1 && s.amplicons <= 2);
            note = "ローカルで予測PCR産物が2本以下（1本のみはなし）";
          }

          if (picked.length > 0) {
            picked.sort((a, b) => {
              if (a.amplicons !== b.amplicons) return a.amplicons - b.amplicons;
              if (a.totalHits !== b.totalHits) return a.totalHits - b.totalHits;
              return (a.pair.pair_penalty ?? 0) - (b.pair.pair_penalty ?? 0);
            });
            filteredCandidates = picked.slice(0, Math.max(1, autoBlastTopN)).map((s) => ({
              pair: s.pair,
              left: s.left,
              right: s.right,
              note,
              ampliconCount: s.amplicons,
              quality: gradeFromAmp(s.amplicons),
            }));
          } else {
            // 予測PCR産物がどのペアでも見つからない場合は、
            // 旧挙動と同様に filteredCandidates を使わず全候補を UI に出す。
            filteredCandidates = undefined;
          }
        }

        return { ...plan, range: [s, e], result: res, filteredCandidates, loading: false };
      } catch (err) {
        const msg = err instanceof Error ? err.message : "エラーが発生しました";
        return { ...plan, error: msg, loading: false };
      }
    };

    let nextIndex = 0;
    const worker = async () => {
      for (; ;) {
        if (batchCancelRef.current) return;
        const idx = nextIndex;
        nextIndex += 1;
        if (idx >= plansToDesign.length) return;
        const plan = plansToDesign[idx];
        const updatedPlan = await runPlan(plan);
        setExonPrimerPlans((prev) => {
          const next = prev.slice();
          next[idx] = updatedPlan;
          return next;
        });
        setBatchStatus((prev) =>
          prev ? { ...prev, done: Math.min(prev.total, prev.done + 1) } : prev,
        );
      }
    };

    await Promise.all(Array.from({ length: parallel }, worker));

    if (batchCancelRef.current) {
      setExonPrimerPlans((prev) =>
        prev.map((p) => (p.loading ? { ...p, loading: false, error: "キャンセルしました" } : p)),
      );
    }

    setBatchStatus((prev) =>
      prev ? { ...prev, cancelled: batchCancelRef.current } : prev,
    );
  };

  const renderBlastTable = (res: BlastResponse | null, title: string) => {
    if (!res) return null;
    if (res.hits.length === 0) return <p className="seq-hint">{title}: ヒットなし</p>;
    const grouped = res.hits.reduce<Record<string, BlastHit[]>>((acc, h) => {
      const src = h.source || "ncbi";
      acc[src] = acc[src] || [];
      acc[src].push(h);
      return acc;
    }, {});
    return (
      <div>
        <h4>{title}</h4>
        {Object.entries(grouped).map(([src, hits]) => (
          <div key={src} className="primer-blast-table">
            <p className="seq-hint">{src} ({hits.length} hits)</p>
            <div className="table-scroll">
              <table className="seq-table">
                <thead>
                  <tr>
                    <th>#</th>
                    <th>ヒット ID</th>
                    <th>% identity</th>
                    <th>長さ</th>
                    <th>E-value</th>
                    <th>範囲</th>
                    <th>source</th>
                  </tr>
                </thead>
                <tbody>
                  {hits.map((hit, idx) => (
                    <tr key={`${hit.sseqid}-${idx}`}>
                      <td>{idx + 1}</td>
                      <td>{hit.sseqid}</td>
                      <td>{hit.pident.toFixed(1)}</td>
                      <td>{hit.length}</td>
                      <td>{hit.evalue.toExponential(2)}</td>
                      <td>
                        {hit.qstart}–{hit.qend}
                      </td>
                      <td>{hit.source ?? "-"}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        ))}
      </div>
    );
  };

  const designPrimersForCds = async () => {
    if (!normalizedSeq) {
      setDesignError("配列が読み込まれていません。");
      return;
    }
    if (parsedCds.length === 0) {
      setDesignError("CDS 情報がありません。CDS 範囲または CDS 配列を確認してください。");
      return;
    }
    setDesignError(null);

    const useSequencingLogic = mode === "sequence" ? true : sequencingMode;

    let plansToDesign: Array<{
      kind?: "exon" | "cds";
      exonIndex: number;
      exonEndIndex?: number;
      subIndex?: number;
      range: [number, number];
      result?: PrimerDesignResponse;
      error?: string;
      loading?: boolean;
      filteredCandidates?: Array<{
        pair: PrimerPair;
        left?: BlastResponse;
        right?: BlastResponse;
        note?: string;
        ampliconCount?: number;
        quality?: "S" | "A" | "B" | "C" | "D";
      }>;
    }> = [];

    const addTiledRanges = (start: number, end: number, cdsIndex: number) => {
      let currentStart = start;
      let subIdx = 1;
      while (currentStart < end) {
        let currentEnd = currentStart + seqProductMax - 1;
        if (currentEnd > end) currentEnd = end;

        if (currentEnd - currentStart + 1 < seqProductMin) {
          currentStart = Math.max(start, currentEnd - seqProductMin + 1);
        }

        plansToDesign.push({
          kind: "cds",
          exonIndex: cdsIndex,
          subIndex: subIdx,
          range: [currentStart, currentEnd],
          loading: true,
        });

        if (currentEnd >= end) break;
        currentStart = currentEnd - seqOverlap + 1;
        subIdx += 1;
      }
    };

    if (useSequencingLogic) {
      // CDS 全体を 600–800bp のタイルでカバーする。
      const n = parsedCds.length;
      let i = 0;

      while (i < n) {
        const [sRaw, eRaw] = parsedCds[i];
        const startWithMargin = Math.max(1, sRaw - seqMargin);
        const endWithMargin = Math.min(normalizedSeq.length, eRaw + seqMargin);
        const spanLen = endWithMargin - startWithMargin + 1;

        if (spanLen > seqProductMax) {
          // この CDS セグメント単体が長すぎるので、セグメントごとにタイル分割
          addTiledRanges(startWithMargin, endWithMargin, i + 1);
          i += 1;
          continue;
        }

        // まず現在の CDS セグメントを含むブロックを開始
        let blockStart = startWithMargin;
        let blockEnd = endWithMargin;
        const blockCdsIndex = i + 1;
        let blockCdsEndIndex = i + 1;
        let j = i + 1;

        // 可能なら次の CDS セグメントも同じブロックにまとめる（産物長 <= max）
        while (j < n) {
          const [nextS, nextE] = parsedCds[j];
          const nextStart = Math.max(1, nextS - seqMargin);
          const nextEnd = Math.min(normalizedSeq.length, nextE + seqMargin);
          const newSpan = nextEnd - blockStart + 1;
          if (newSpan <= seqProductMax) {
            blockEnd = nextEnd;
            blockCdsEndIndex = j + 1;
            j += 1;
          } else {
            break;
          }
        }

        const blockLen = blockEnd - blockStart + 1;
        if (blockLen <= seqProductMax) {
          plansToDesign.push({
            kind: "cds",
            exonIndex: blockCdsIndex,
            exonEndIndex: blockCdsEndIndex,
            range: [blockStart, blockEnd],
            loading: true,
          });
        } else {
          addTiledRanges(blockStart, blockEnd, blockCdsIndex);
        }

        i = j;
      }
    } else {
      // Sequencing Mode でない場合は CDS セグメントごとに 1 本ずつ設計
      plansToDesign = parsedCds.map((r, idx) => ({
        kind: "cds",
        exonIndex: idx + 1,
        range: r,
        loading: true,
      }));
    }

    setExonPrimerPlans(plansToDesign);
    batchCancelRef.current = false;
    setBatchStatus(
      plansToDesign.length
        ? { total: plansToDesign.length, done: 0, cancelled: false }
        : null,
    );

    const hw = typeof navigator !== "undefined" ? navigator.hardwareConcurrency || 12 : 12;
    const baseParallel = Math.max(1, Math.min(12, Math.floor(hw / 2) || 1));
    const parallel = Math.max(
      1,
      Math.min(
        plansToDesign.length,
        autoBlastExons ? Math.min(4, baseParallel) : baseParallel,
      ),
    );

    const runPlan = async (
      plan: (typeof plansToDesign)[number],
    ): Promise<(typeof plansToDesign)[number]> => {
      let [s, e] = plan.range;
      try {
        const len = Math.max(1, e - s + 1);
        let productRange = "";

        if (useSequencingLogic) {
          productRange = `${seqProductMin}-${seqProductMax}`;
        } else {
          const lower = Math.max(80, len);
          const upper = Math.min(normalizedSeq.length, lower + 150);
          productRange = `${lower}-${upper}`;
        }

        let res = await bioapiClient.designPrimers({
          sequence: normalizedSeq,
          num_return: Math.min(numReturn, 5),
          product_size_range: productRange,
          target_start: s,
          target_length: len,
          opt_tm: pcrOptTm,
          min_tm: pcrMinTm,
          max_tm: pcrMaxTm,
          primer_min_size: pcrMinSize,
          primer_opt_size: pcrOptSize,
          primer_max_size: pcrMaxSize,
          primer_min_gc: pcrMinGc,
          primer_max_gc: pcrMaxGc,
          primer_salt_monovalent: pcrSaltMonovalent,
          primer_dna_conc: pcrDnaConc,
        });

        if (useSequencingLogic && !res.candidates.length) {
          const rangesForKind = plan.kind === "cds" ? parsedCds : parsedExons;
          const span = rangesForKind[plan.exonIndex - 1];
          if (span) {
            const [spanStart, spanEnd] = span;
            const allowStart = Math.max(1, spanStart - seqMargin);
            const allowEnd = Math.min(normalizedSeq.length, spanEnd + seqMargin);
            const maxWindow = seqProductMax;
            const minWindow = seqProductMin;

            const lengthsToTry: number[] = [];
            const initialLen = Math.min(len, maxWindow);
            const minTargetLen = Math.max(
              Math.min(len, minWindow),
              Math.floor(minWindow * 0.6),
            );
            lengthsToTry.push(initialLen);
            for (let t = initialLen - 50; t >= minTargetLen; t -= 50) {
              if (!lengthsToTry.includes(t)) lengthsToTry.push(t);
            }

            let best: PrimerDesignResponse | null = null;
            let bestStart = s;
            let bestLen = len;
            let trialCount = 0;
            const spanLength = allowEnd - allowStart + 1;
            const step = 30;
            const approxPerLen = Math.max(1, Math.ceil(spanLength / step) * 2);
            const maxTrials = Math.min(800, approxPerLen * Math.max(1, lengthsToTry.length));

            outer: for (const tLen of lengthsToTry) {
              const maxStart = allowEnd - tLen + 1;
              const minStart = allowStart;
              if (maxStart < minStart) {
                continue;
              }
              const centerStart = Math.min(Math.max(s, minStart), maxStart);

              for (let offset = 0; ; offset += step) {
                const candidateStarts =
                  offset === 0 ? [centerStart] : [centerStart - offset, centerStart + offset];
                for (const newStart of candidateStarts) {
                  if (newStart < minStart || newStart > maxStart) continue;
                  const newEnd = newStart + tLen - 1;
                  if (newEnd > allowEnd) continue;

                  trialCount += 1;
                  if (trialCount > maxTrials) {
                    break outer;
                  }

                  const trial = await bioapiClient.designPrimers({
                    sequence: normalizedSeq,
                    num_return: Math.min(numReturn, 5),
                    product_size_range: productRange,
                    target_start: newStart,
                    target_length: tLen,
                    opt_tm: pcrOptTm,
                    min_tm: pcrMinTm,
                    max_tm: pcrMaxTm,
                    primer_min_size: pcrMinSize,
                    primer_opt_size: pcrOptSize,
                    primer_max_size: pcrMaxSize,
                    primer_min_gc: pcrMinGc,
                    primer_max_gc: pcrMaxGc,
                    primer_salt_monovalent: pcrSaltMonovalent,
                    primer_dna_conc: pcrDnaConc,
                  });
                  if (trial.candidates.length) {
                    best = trial;
                    bestStart = newStart;
                    bestLen = tLen;
                    break outer;
                  }
                }
                if (centerStart - offset <= minStart && centerStart + offset >= maxStart) {
                  break;
                }
              }
            }

            if (best) {
              res = best;
              s = bestStart;
              e = bestStart + bestLen - 1;
            }
          }
        } else if (!useSequencingLogic && mode === "exon" && !res.candidates.length) {
          const relaxed = await bioapiClient.designPrimers({
            sequence: normalizedSeq,
            num_return: Math.min(numReturn, 5),
            product_size_range: null,
            target_start: s,
            target_length: len,
            opt_tm: pcrOptTm,
            min_tm: pcrMinTm,
            max_tm: pcrMaxTm,
            primer_min_size: pcrMinSize,
            primer_opt_size: pcrOptSize,
            primer_max_size: pcrMaxSize,
            primer_min_gc: pcrMinGc,
            primer_max_gc: pcrMaxGc,
            primer_salt_monovalent: pcrSaltMonovalent,
            primer_dna_conc: pcrDnaConc,
          });
          if (relaxed.candidates.length) {
            res = relaxed;
          }
        }

        let filteredCandidates:
          | Array<{
            pair: PrimerPair;
            left: BlastResponse;
            right: BlastResponse;
            note?: string;
            ampliconCount?: number;
          }>
          | undefined;

        if (autoBlastExons) {
          const scored: Array<{
            pair: PrimerPair;
            left: BlastResponse;
            right: BlastResponse;
            totalHits: number;
            amplicons: number;
          }> = [];
          const candidatesToEval = res.candidates.slice(
            0,
            Math.max(autoBlastTopN * 2, autoBlastTopN, 3),
          );
          if (blastUseLocal && effectiveLocalDbs.length > 0) {
            const idxBySeq = new Map<string, number>();
            const uniq: string[] = [];
            const normalize = (s: string) => s.replace(/\s+/g, "").toUpperCase();
            const indexOf = (s: string) => {
              const key = normalize(s);
              const existing = idxBySeq.get(key);
              if (existing != null) return existing;
              const idx = uniq.length;
              idxBySeq.set(key, idx);
              uniq.push(key);
              return idx;
            };

            candidatesToEval.forEach((p) => {
              indexOf(p.left_sequence);
              indexOf(p.right_sequence);
            });

            if (uniq.length > 0) {
              const batch = await bioapiClient.runBlastBatchLocal({
                sequences: uniq,
                dbs: effectiveLocalDbs,
                local_mode: localMode,
                task: blastTask,
                evalue: blastEvalue,
                max_target_seqs: blastMaxHits,
                max_hsps: blastMaxHsps ?? undefined,
                num_threads: blastNumThreads ?? undefined,
              });
              if (!batch.results || batch.results.length !== uniq.length) {
                throw new Error("run_batch_local の結果件数が期待値と一致しませんでした。");
              }

              const perSeq = batch.results.map((r) => ({
                ...r,
                hits: relabelLocalBlastHits(r.hits ?? [], effectiveLocalDbs, localDbOptions),
              }));

              for (const pair of candidatesToEval) {
                const left = perSeq[indexOf(pair.left_sequence)];
                const right = perSeq[indexOf(pair.right_sequence)];
                const totalHits = countLocalHits(left) + countLocalHits(right);
                const { amplicons } = computePrimerAmplicons(left, right);
                scored.push({ pair, left, right, totalHits, amplicons: amplicons.length });
              }
            }
          } else {
            for (const pair of candidatesToEval) {
              const [left, right] = await Promise.all([
                runBlastForSequence(pair.left_sequence),
                runBlastForSequence(pair.right_sequence),
              ]);
              const totalHits = countLocalHits(left) + countLocalHits(right);
              const { amplicons } = computePrimerAmplicons(left, right);
              scored.push({ pair, left, right, totalHits, amplicons: amplicons.length });
            }
          }

          const strict = scored.filter((s) => s.amplicons === 1);
          let picked = strict;
          let note: string | undefined = "ローカルで予測PCR産物が1本のみ";

          if (picked.length === 0) {
            picked = scored.filter((s) => s.amplicons >= 1 && s.amplicons <= 2);
            note = "ローカルで予測PCR産物が2本以下（1本のみはなし）";
          }

          if (picked.length > 0) {
            picked.sort((a, b) => {
              if (a.amplicons !== b.amplicons) return a.amplicons - b.amplicons;
              if (a.totalHits !== b.totalHits) return a.totalHits - b.totalHits;
              return (a.pair.pair_penalty ?? 0) - (b.pair.pair_penalty ?? 0);
            });
            filteredCandidates = picked.slice(0, Math.max(1, autoBlastTopN)).map((s) => ({
              pair: s.pair,
              left: s.left,
              right: s.right,
              note,
              ampliconCount: s.amplicons,
            }));
          } else {
            filteredCandidates = undefined;
          }
        }

        return { ...plan, range: [s, e], result: res, filteredCandidates, loading: false };
      } catch (err) {
        const msg = err instanceof Error ? err.message : "エラーが発生しました";
        return { ...plan, error: msg, loading: false };
      }
    };

    let nextIndex = 0;
    const worker = async () => {
      for (; ;) {
        if (batchCancelRef.current) return;
        const idx = nextIndex;
        nextIndex += 1;
        if (idx >= plansToDesign.length) return;
        const plan = plansToDesign[idx];
        const updatedPlan = await runPlan(plan);
        setExonPrimerPlans((prev) => {
          const next = prev.slice();
          next[idx] = updatedPlan;
          return next;
        });
        setBatchStatus((prev) =>
          prev ? { ...prev, done: Math.min(prev.total, prev.done + 1) } : prev,
        );
      }
    };

    await Promise.all(Array.from({ length: parallel }, worker));

    if (batchCancelRef.current) {
      setExonPrimerPlans((prev) =>
        prev.map((p) => (p.loading ? { ...p, loading: false, error: "キャンセルしました" } : p)),
      );
    }

    setBatchStatus((prev) =>
      prev ? { ...prev, cancelled: batchCancelRef.current } : prev,
    );
  };

  const summarizeBlastHitsPerDb = (
    left: BlastResponse | undefined,
    right: BlastResponse | undefined,
  ): Array<{ db: string; left: number; right: number }> => {
    const map = new Map<string, { left: number; right: number }>();
    const addHits = (resp: BlastResponse | undefined, side: "left" | "right") => {
      if (!resp) return;
      resp.hits.forEach((h) => {
        const raw = h.source || (side === "left" ? "left" : "right");
        const dbLabel = raw.startsWith("local:") ? raw.replace(/^local:/, "") : raw;
        const rec = map.get(dbLabel) || { left: 0, right: 0 };
        rec[side] += 1;
        map.set(dbLabel, rec);
      });
    };
    addHits(left, "left");
    addHits(right, "right");

    // ローカル DB については、ヒット 0 のものも明示的に表示
    effectiveLocalDbs.forEach((path) => {
      const label = labelForDb(path);
      if (!map.has(label)) {
        map.set(label, { left: 0, right: 0 });
      }
    });

    return Array.from(map.entries())
      .map(([db, counts]) => ({ db, ...counts }))
      .sort((a, b) => a.db.localeCompare(b.db));
  };

  const buildMarkdownReport = (): string => {
    if (!normalizedSeq || exonPrimerPlans.length === 0) return "";

    const dt = new Date();
    const headerLines: string[] = [];
    headerLines.push("# Exon/CDS プライマー設計レポート");
    headerLines.push("");
    headerLines.push(`- 作成時刻: ${dt.toLocaleString()}`);
    if (geneId) headerLines.push(`- Gene ID: \`${geneId}\``);
    if (species) headerLines.push(`- 種: \`${species}\``);
    headerLines.push(
      `- モード: ${mode === "sequence" || sequencingMode
        ? "シーケンスプライマー（600–800bp タイル設計）"
        : "CDS/エキソン増幅"
      }`,
    );
    if (mode === "sequence" || sequencingMode) {
      headerLines.push(
        `- シーケンス産物長: \`${seqProductMin}-${seqProductMax} bp\` ／ オーバーラップ: \`${seqOverlap} bp\` ／ マージン: \`${seqMargin} bp\``,
      );
    }
    if (exonRanges) headerLines.push(`- エキソン範囲（ゲノム座標）: \`${exonRanges}\``);
    if (cdsRange) headerLines.push(`- CDS 範囲（ゲノム座標）: \`${cdsRange}\``);
    headerLines.push("");

    const lines: string[] = [...headerLines];
    lines.push("## サマリ");
    lines.push("");
    lines.push(
      "| ターゲット | ゲノム範囲 | 産物長 (bp) | 品質ランク (S–D) | Primer3 penalty |",
    );
    lines.push("| --- | --- | ---: | :--: | ---: |");

    exonPrimerPlans.forEach((plan) => {
      const top =
        plan.filteredCandidates?.[0] ??
        (plan.result && plan.result.candidates[0]
          ? { pair: plan.result.candidates[0] }
          : null);
      if (!top) return;
      const cand = (top as any).pair as PrimerPair;
      const quality = (top as any).quality as
        | "S"
        | "A"
        | "B"
        | "C"
        | "D"
        | undefined;
      const label =
        (plan.kind ?? "exon") === "cds"
          ? `CDS${plan.exonIndex}`
          : `Exon${plan.exonIndex}`;
      const rangeLabel = `${plan.range[0]}–${plan.range[1]}`;
      const prod = cand.product_size ?? "";
      const penalty =
        typeof cand.pair_penalty === "number"
          ? cand.pair_penalty.toFixed(2)
          : "";
      lines.push(
        `| ${label} | ${rangeLabel} | ${prod} | ${quality ?? ""} | ${penalty} |`,
      );
    });

    lines.push("");
    lines.push("## ターゲット別詳細");
    lines.push("");

    exonPrimerPlans.forEach((plan) => {
      const titleLabel =
        (plan.kind ?? "exon") === "cds"
          ? `CDS ${plan.exonIndex}`
          : `Exon ${plan.exonIndex}`;
      lines.push(
        `### ${titleLabel} (${plan.range[0]}–${plan.range[1]} bp)${plan.subIndex ? ` サブ領域 ${plan.subIndex}` : ""
        }`,
      );
      lines.push("");
      if (!plan.result || plan.result.candidates.length === 0) {
        lines.push("_候補が見つかりませんでした。_");
        lines.push("");
        return;
      }

      const rows =
        plan.filteredCandidates && plan.filteredCandidates.length > 0
          ? plan.filteredCandidates
          : plan.result.candidates.map((c) => ({ pair: c }));

      lines.push(
        "| # | 左プライマー | 右プライマー | 産物長 (bp) | 品質ランク | Primer3 penalty | 予測 PCR 産物数 (ローカル) | 備考 |",
      );
      lines.push(
        "| ---: | --- | --- | ---: | :--: | ---: | ---: | --- |",
      );

      rows.forEach((entry, idx) => {
        const cand = (entry as any).pair as PrimerPair;
        const quality = (entry as any).quality as
          | "S"
          | "A"
          | "B"
          | "C"
          | "D"
          | undefined;
        const ampCount = (entry as any).ampliconCount as
          | number
          | undefined;
        const note = (entry as any).note as string | undefined;
        const penalty =
          typeof cand.pair_penalty === "number"
            ? cand.pair_penalty.toFixed(2)
            : "";

        lines.push(
          `| ${idx + 1} | \`${cand.left_sequence}\` | \`${cand.right_sequence}\` | ${cand.product_size ?? ""
          } | ${quality ?? ""} | ${penalty} | ${typeof ampCount === "number" ? ampCount : ""
          } | ${note ?? ""} |`,
        );

        const left = (entry as any).left as BlastResponse | undefined;
        const right = (entry as any).right as BlastResponse | undefined;
        if (left || right) {
          const perDb = summarizeBlastHitsPerDb(left, right);
          if (perDb.length) {
            lines.push("");
            lines.push("DB 別 BLAST ヒット数:");
            perDb.forEach((d) => {
              lines.push(
                `- ${d.db}: Forward ${d.left} 件 / Reverse ${d.right} 件`,
              );
            });
          }
        }
      });

      lines.push("");
    });

    if (currentPatternPrimers.length > 0) {
      lines.push("## Primer 逆引き用ペア一覧");
      lines.push("");
      lines.push(
        "以下のブロックを Primer 逆引きタブの「プライマー配列（2 行ずつ 1 ペア）」に貼り付けてください。",
      );
      lines.push("");
      lines.push("```");
      currentPatternPrimers.forEach((p) => {
        lines.push(p.left);
        lines.push(p.right);
      });
      lines.push("```");
      lines.push("");
    }

    return lines.join("\n");
  };

  const autoMapFromSequences = () => {
    if (!normalizedSeq) {
      setStructureError("先にゲノム配列を貼り付けてください。");
      return;
    }
    setStructureError(null);
    setInfoMessage("");
    const genome = normalizedSeq;

    const exonSeqs = exonSequencesInput
      .split(/\r?\n/)
      .map((s) => s.trim())
      .filter(Boolean)
      .map((s) => s.replace(/\s+/g, "").toUpperCase());

    const cdsSeq = cdsSequenceInput.replace(/\s+/g, "").toUpperCase();

    const ranges: string[] = [];
    let cdsMapped = false;
    let mappedAny = false;
    for (const ex of exonSeqs) {
      if (!ex) continue;
      let pos = genome.indexOf(ex);
      let length = ex.length;
      if (pos < 0) {
        const rc = revComp(ex);
        pos = genome.indexOf(rc);
      }
      if (pos < 0) {
        setStructureError(`エキソン配列がゲノム中に見つかりませんでした: ${ex.slice(0, 20)}...`);
        return;
      }
      const start = pos + 1;
      const end = pos + length;
      ranges.push(`${start}-${end}`);
      mappedAny = true;
    }
    if (ranges.length > 0) {
      // ソートして連続表示
      const parts = ranges
        .map((r) => r.split("-").map((n) => Number(n)))
        .filter((p) => p.length === 2 && Number.isFinite(p[0]) && Number.isFinite(p[1]))
        .map(([a, b]) => [Math.min(a, b), Math.max(a, b)] as [number, number])
        .sort((a, b) => a[0] - b[0]);
      setExonRanges(parts.map(([s, e]) => `${s}-${e}`).join(","));
      const minStart = Math.min(...parts.map((p) => p[0]));
      const maxEnd = Math.max(...parts.map((p) => p[1]));
      setSelectionRange({ start: minStart, end: maxEnd });
      setSubStart(minStart);
      setSubEnd(maxEnd);
      mappedAny = true;
    }

    if (cdsSeq) {
      let pos = genome.indexOf(cdsSeq);
      let len = cdsSeq.length;
      if (pos < 0) {
        const rc = revComp(cdsSeq);
        pos = genome.indexOf(rc);
      }
      if (pos >= 0) {
        const start = pos + 1;
        const end = pos + len;
        setCdsRange(`${start}-${end}`);
        cdsMapped = true;
        mappedAny = true;
      } else if (ranges.length > 0) {
        // cDNA/CDS がスプライス済みならエキソン範囲をそのまま CDS として扱う
        setCdsRange(ranges.join(","));
        cdsMapped = true;
      }
    } else if (ranges.length > 0) {
      // fallback to exon span
      const minStart = Math.min(...ranges.map((r) => Number(r.split("-")[0])));
      const maxEnd = Math.max(...ranges.map((r) => Number(r.split("-")[1])));
      setCdsRange(`${minStart}-${maxEnd}`);
    }

    if (mappedAny) {
      setInfoMessage("ゲノム配列上にエキソン/CDS をマッピングしました。");
    } else {
      setStructureError("貼り付けた配列がゲノム中に見つかりませんでした。");
    }
  };

  const parseFastaBlocks = (text: string) => {
    const entries: { header: string; seq: string }[] = [];
    let currentHeader: string | null = null;
    let currentSeq: string[] = [];
    const lines = text.split(/\r?\n/);
    for (const raw of lines) {
      const line = raw.trim();
      if (!line) continue;
      if (line.startsWith(">")) {
        if (currentHeader) {
          entries.push({ header: currentHeader, seq: currentSeq.join("") });
        }
        currentHeader = line.slice(1);
        currentSeq = [];
      } else {
        currentSeq.push(line.replace(/\s+/g, "").toUpperCase());
      }
    }
    if (currentHeader) {
      entries.push({ header: currentHeader, seq: currentSeq.join("") });
    }
    return entries;
  };

  const computeExonRangesOnSeq = (seq: string, exonSeqs: string[]) => {
    const ranges: string[] = [];
    for (const ex of exonSeqs) {
      if (!ex) continue;
      let pos = seq.indexOf(ex);
      let length = ex.length;
      if (pos < 0) {
        const rc = revComp(ex);
        pos = seq.indexOf(rc);
      }
      if (pos < 0) return null;
      const start = pos + 1;
      const end = pos + length;
      ranges.push(`${start}-${end}`);
    }
    return ranges.join(",");
  };

  const loadFromFasta = () => {
    const blocks = parseFastaBlocks(fastaInput);
    if (!blocks.length) {
      setStructureError("FASTA を貼り付けてください。");
      return;
    }
    const header0 = blocks[0].header;
    const inferredGene = header0.split(/\s+/)[0];
    const findEntry = (kw: RegExp) =>
      blocks.find((b) => kw.test(b.header.toLowerCase()));

    const chooseGenome = () => {
      const byKeyword =
        blocks.find((b) =>
          /primary_assembly|assembly|genome|genomic|gdna|dna[:\s]|chr|lg|cm\d/i.test(
            b.header.toLowerCase(),
          ),
        ) || null;
      if (byKeyword) return byKeyword;
      // キーワードが無ければ最長配列をゲノム候補とする
      return blocks.reduce((acc, cur) => (cur.seq.length > acc.seq.length ? cur : acc), blocks[0]);
    };

    const genome = chooseGenome();
    const cDNA = findEntry(/cdna/);
    const cds = findEntry(/cds/);
    const exonEntries = blocks.filter(
      (b) => /exon/.test(b.header.toLowerCase()) || /-e\d+/i.test(b.header),
    );

    // ベースのゲノム配列は textarea に貼られたものを優先
    let genomeSeq = normalizedSeq;
    if (!genomeSeq && genome?.seq) {
      genomeSeq = genome.seq.replace(/\s+/g, "");
      setSequence(genome.seq);
      setGenomeHeader(genome.header || "genome_sequence");
    }
    if (!genomeSeq) {
      setStructureError("先にゲノムDNA配列を貼り付けてください（または FASTA に genome/assembly を含めてください）。");
      return;
    }
    const exonSeqs = exonEntries.map((e) => e.seq.replace(/\s+/g, ""));

    setGeneId((prev) => (prev ? prev : inferredGene || ""));
    setCdsSequenceInput(cds?.seq || cDNA?.seq || "");
    setExonSequencesInput(exonSeqs.join("\n"));

    const exonRangesStr =
      exonSeqs.length > 0 ? computeExonRangesOnSeq(genomeSeq, exonSeqs) : null;
    if (exonRangesStr) {
      setExonRanges(exonRangesStr);
    } else if (exonSeqs.length > 0) {
      setStructureError("エキソン配列をゲノム上で位置決めできませんでした。順番・配列を確認してください。");
      return;
    }
    const cdsSeq = (cds?.seq || cDNA?.seq || "").replace(/\s+/g, "");
    const cdsLen = cdsSeq.length;
    if (exonRangesStr) {
      const parts = exonRangesStr
        .split(",")
        .map((r) => r.split("-").map((n) => Number(n)))
        .filter((p) => p.length === 2 && Number.isFinite(p[0]) && Number.isFinite(p[1]));
      const minStart = Math.min(...parts.map((p) => Math.min(...p)));
      const maxEnd = Math.max(...parts.map((p) => Math.max(...p)));
      // CDS がスプライス済みでゲノム上で連続しないため、エキソン範囲を CDS として保持
      setCdsRange(exonRangesStr);
      setSelectionRange({ start: minStart, end: maxEnd });
      setSubStart(minStart);
      setSubEnd(maxEnd);
    } else {
      if (cdsSeq) {
        let pos = genomeSeq.indexOf(cdsSeq);
        if (pos < 0) {
          const rc = revComp(cdsSeq);
          pos = genomeSeq.indexOf(rc);
        }
        if (pos >= 0) {
          const start = pos + 1;
          const end = pos + cdsSeq.length;
          setCdsRange(`${start}-${end}`);
          setSelectionRange({ start, end });
          setSubStart(start);
          setSubEnd(end);
        } else {
          setCdsRange(`1-${genomeSeq.length}`);
          setSelectionRange(null);
        }
      } else {
        setCdsRange(`1-${genomeSeq.length}`);
        setSelectionRange(null);
      }
    }
    setInfoMessage(
      genomeSeq
        ? "ゲノムDNAを基準に cDNA / CDS / Exon を自動設定しました。"
        : "ゲノムDNAが空です。先にゲノム配列を貼り付けてください。",
    );
    setStructureError(null);
  };

  const loadSample = () => {
    setCdsSequenceInput(sampleCds);
    setExonSequencesInput(sampleExons);
    setInfoMessage("サンプルの CDS / エキソン配列をセットしました。ゲノム配列を貼り付けてから自動算出を押してください。");
  };

  return (
    <div className="exon-grid">
      <div className="exon-controls">
        <h3>
          {mode === "sequence"
            ? "シーケンスプライマー設計（エキソンカバー）"
            : "CDS/エキソン ハイライト & プライマー設計"}
        </h3>

        {/* --- Gene Loader Section --- */}
        <div className="primer-row" style={{ alignItems: "flex-end", marginBottom: "0.8rem", padding: "0.5rem", border: "1px solid #ddd", borderRadius: "4px" }}>
          <label className="seq-label" style={{ width: "auto" }}>
            Target DB:
            <select
              className="seq-input"
              style={{ minWidth: "140px" }}
              value={loadGeneDb}
              onChange={(e) => setLoadGeneDb(e.target.value)}
              disabled={loadGeneLoading}
            >
              <option value="auto">Auto Detect (All)</option>
              {localDbOptions.map((o) => (
                <option key={o.value} value={o.value}>
                  {o.label}
                </option>
              ))}
              <option value="custom">Custom...</option>
            </select>
          </label>

          {loadGeneDb === "custom" && (
            <label className="seq-label" style={{ width: "auto" }}>
              Custom Path:
              <input
                className="seq-input"
                value={loadGeneCustomDb}
                onChange={(e) => setLoadGeneCustomDb(e.target.value)}
                placeholder="/path/to/db"
              />
            </label>
          )}

          <label className="seq-label">
            Gene Name (e.g. GENE...):
            <input
              className="seq-input"
              value={loadGeneName}
              onChange={(e) => setLoadGeneName(e.target.value)}
              placeholder="Gene ID / Name"
              onKeyDown={(e) => {
                if (e.key === "Enter") handleLoadGene();
              }}
              disabled={loadGeneLoading}
            />
          </label>

          <label className="seq-label" style={{ width: "100px" }}>
            Margin (bp):
            <input
              className="seq-input"
              type="number"
              value={loadGeneMargin}
              onChange={(e) => {
                const v = e.target.value;
                setLoadGeneMargin(v === "" ? "" : Number(v));
              }}
              placeholder="Auto"
              disabled={loadGeneLoading}
            />
          </label>

          <button
            type="button"
            className="seq-button"
            onClick={handleLoadGene}
            disabled={loadGeneLoading || !loadGeneName}
          >
            {loadGeneLoading ? "Loading..." : "Load Gene"}
          </button>
        </div>
        {/* ------------------------- */}

        <div className="dual-input-row">
          <div>
            <label className="seq-label" style={{ flex: 1 }}>
              FASTA 形式で一括入力 (DNA/CDS/Exon):
              <textarea
                className="seq-textarea"
                rows={5}
                value={fastaInput}
                onChange={(e) => setFastaInput(e.target.value)}
                placeholder=">id cdna&#10;ATGC...&#10;>id cds&#10;ATGC...&#10;>id exon1&#10;..."
              />
            </label>
            <div className="primer-controls" style={{ minWidth: "180px" }}>
              <button type="button" className="seq-button secondary" onClick={loadFromFasta}>
                FASTA からセット
              </button>
              <p className="seq-hint">
                cdna / cds / exon を含むヘッダを自動検出し、CDS 範囲とエキソンを配置します。
              </p>
            </div>
          </div>

          <div>
            <label className="seq-label">
              ゲノムDNA配列:
              <textarea
                className="seq-textarea"
                rows={8}
                value={sequence}
                onChange={(e) => {
                  const val = e.target.value;
                  const headerLine = val
                    .split(/\r?\n/)
                    .find((l) => l.trim().startsWith(">"));
                  if (headerLine) {
                    setGenomeHeader(headerLine.replace(/^>/, "").trim() || "genome_sequence");
                  }
                  setSequence(val);
                }}
                placeholder="ゲノムDNA配列を貼り付け"
              />
            </label>
          </div>
        </div>
        {mode === "exon" && (
          <>
            <div className="primer-row">
              <label className="seq-label">
                Feature ID:
                <input
                  className="seq-input"
                  type="text"
                  value={geneId}
                  onChange={(e) => setGeneId(e.target.value)}
                  placeholder="例: feature_001"
                />
              </label>
              <label className="seq-label">
                species:
                <input
                  className="seq-input"
                  type="text"
                  value={species}
                  onChange={(e) => setSpecies(e.target.value)}
                  placeholder="e.g. arabidopsis_thaliana"
                />
              </label>
            </div>
            <div className="primer-row">
              <button
                type="button"
                className="seq-button secondary"
                onClick={fetchStructure}
                disabled={structureLoading}
              >
                {structureLoading ? "取得中..." : "構造を取得"}
              </button>
              <button
                type="button"
                className="seq-button secondary"
                onClick={guessGeneIdByBlast}
                disabled={structureLoading}
              >
                {structureLoading ? "取得中..." : "BLAST で Gene 候補を推定"}
              </button>
            </div>
            {structureLoading ? (
              <JobProgressCard title="ローカル BLAST" jobId={blastJobId} job={blastJobInfo} />
            ) : null}
          </>
        )}
        <div className="primer-row">
          <label className="seq-label">
            CDS 配列（任意、spliced CDS）:
            <textarea
              className="seq-textarea"
              rows={3}
              value={cdsSequenceInput}
              onChange={(e) => setCdsSequenceInput(e.target.value)}
              placeholder="spliced CDS 配列を貼り付け"
            />
          </label>
          <label className="seq-label">
            エキソン配列（任意、1行1エキソン）:
            <textarea
              className="seq-textarea"
              rows={3}
              value={exonSequencesInput}
              onChange={(e) => setExonSequencesInput(e.target.value)}
              placeholder="1行1エキソンで貼り付け"
            />
          </label>
        </div>
        {mode === "exon" && (
          <div className="primer-row">
            <label className="seq-label" style={{ fontWeight: "bold", color: "#0d47a1" }}>
              <input
                type="checkbox"
                checked={sequencingMode}
                onChange={(e) => setSequencingMode(e.target.checked)}
                style={{ marginRight: "0.5rem" }}
              />
              Sequencing Mode (Tiling & Margins)
            </label>
          </div>
        )}
        {(mode === "sequence" || sequencingMode) && (
          <div
            className="primer-row"
            style={{ background: "#e3f2fd", padding: "0.5rem", borderRadius: "4px", marginBottom: "1rem" }}
          >
            <label className="seq-label">
              Product Size (bp):
              <div style={{ display: "flex", gap: "0.5rem", alignItems: "center" }}>
                <input
                  className="seq-input"
                  type="number"
                  value={seqProductMin}
                  onChange={(e) => setSeqProductMin(Number(e.target.value))}
                  style={{ width: "80px" }}
                />
                -
                <input
                  className="seq-input"
                  type="number"
                  value={seqProductMax}
                  onChange={(e) => setSeqProductMax(Number(e.target.value))}
                  style={{ width: "80px" }}
                />
              </div>
            </label>
            <label className="seq-label">
              Overlap (bp):
              <input
                className="seq-input"
                type="number"
                value={seqOverlap}
                onChange={(e) => setSeqOverlap(Number(e.target.value))}
                style={{ width: "80px" }}
              />
            </label>
            <label className="seq-label">
              Margin (bp):
              <input
                className="seq-input"
                type="number"
                value={seqMargin}
                onChange={(e) => setSeqMargin(Number(e.target.value))}
                style={{ width: "80px" }}
              />
            </label>
          </div>
        )}
        <div className="primer-row">
          <button type="button" className="seq-button secondary" onClick={autoMapFromSequences}>
            貼り付けた配列から座標を自動算出
          </button>
          <button type="button" className="seq-button secondary" onClick={loadSample}>
            サンプル配列をセット
          </button>
          <button type="button" className="seq-button" onClick={designPrimersForExons}>
            {mode === "sequence" ? "シーケンスプライマーを自動設計" : "各エキソンでプライマー自動設計"}
          </button>
          <button type="button" className="seq-button" onClick={designPrimersForCds}>
            {mode === "sequence" ? "CDS を自動タイル設計" : "CDS でプライマー自動設計"}
          </button>
        </div>
        <div className="primer-row">
          <button
            type="button"
            className="seq-button secondary"
            onClick={() => setShowPcrOptions((v) => !v)}
          >
            {showPcrOptions ? "PCR条件を隠す" : "PCR条件（Tmなど）を表示"}
          </button>
        </div>
        {showPcrOptions && (
          <div
            className="primer-row"
            style={{
              background: "#f9fafb",
              border: "1px dashed #9ca3af",
              borderRadius: "4px",
              padding: "0.6rem",
              flexWrap: "wrap",
              gap: "0.75rem",
            }}
          >
            <label className="seq-label">
              最適 Tm (℃)
              <input
                className="seq-input"
                type="number"
                value={pcrOptTm}
                onChange={(e) => setPcrOptTm(Number(e.target.value) || 60)}
                style={{ width: "90px" }}
              />
            </label>
            <label className="seq-label">
              最小 Tm (℃)
              <input
                className="seq-input"
                type="number"
                value={pcrMinTm}
                onChange={(e) => setPcrMinTm(Number(e.target.value) || 57)}
                style={{ width: "90px" }}
              />
            </label>
            <label className="seq-label">
              最大 Tm (℃)
              <input
                className="seq-input"
                type="number"
                value={pcrMaxTm}
                onChange={(e) => setPcrMaxTm(Number(e.target.value) || 63)}
                style={{ width: "90px" }}
              />
            </label>
            <label className="seq-label">
              プライマー長 min/opt/max
              <div style={{ display: "flex", gap: "0.3rem" }}>
                <input
                  className="seq-input"
                  type="number"
                  value={pcrMinSize}
                  onChange={(e) => setPcrMinSize(Number(e.target.value) || 18)}
                  style={{ width: "70px" }}
                />
                <input
                  className="seq-input"
                  type="number"
                  value={pcrOptSize}
                  onChange={(e) => setPcrOptSize(Number(e.target.value) || 20)}
                  style={{ width: "70px" }}
                />
                <input
                  className="seq-input"
                  type="number"
                  value={pcrMaxSize}
                  onChange={(e) => setPcrMaxSize(Number(e.target.value) || 27)}
                  style={{ width: "70px" }}
                />
              </div>
            </label>
            <label className="seq-label">
              GC% min/max
              <div style={{ display: "flex", gap: "0.3rem" }}>
                <input
                  className="seq-input"
                  type="number"
                  value={pcrMinGc}
                  onChange={(e) => setPcrMinGc(Number(e.target.value) || 20)}
                  style={{ width: "70px" }}
                />
                <input
                  className="seq-input"
                  type="number"
                  value={pcrMaxGc}
                  onChange={(e) => setPcrMaxGc(Number(e.target.value) || 80)}
                  style={{ width: "70px" }}
                />
              </div>
            </label>
            <label className="seq-label">
              Salt (mM)
              <input
                className="seq-input"
                type="number"
                value={pcrSaltMonovalent}
                onChange={(e) => setPcrSaltMonovalent(Number(e.target.value) || 50)}
                style={{ width: "90px" }}
              />
            </label>
            <label className="seq-label">
              DNA 濃度 (nM)
              <input
                className="seq-input"
                type="number"
                value={pcrDnaConc}
                onChange={(e) => setPcrDnaConc(Number(e.target.value) || 50)}
                style={{ width: "90px" }}
              />
            </label>
            <span className="seq-hint">
              Primer3 に渡す Tm / 長さ / GC% / Salt / DNA 濃度です。未指定に戻したい場合は 60 / 57 / 63 と 18/20/27, 20–80%, 50/50 を推奨します。
            </span>
          </div>
        )}
        <div className="primer-row">
          <label className="seq-label" style={{ maxWidth: "320px" }}>
            自動BLASTで特異性フィルタ
            <div style={{ display: "flex", alignItems: "center", gap: "0.4rem" }}>
              <input
                type="checkbox"
                checked={autoBlastExons}
                onChange={(e) => setAutoBlastExons(e.target.checked)}
              />
              <span style={{ fontSize: "0.9rem", color: "#4b5563" }}>上位</span>
              <input
                type="number"
                className="seq-input"
                style={{ width: "70px" }}
                min={1}
                max={5}
                value={autoBlastTopN}
                onChange={(e) => setAutoBlastTopN(Math.max(1, Math.min(5, Number(e.target.value) || 1)))}
                disabled={!autoBlastExons}
              />
              <span style={{ fontSize: "0.9rem", color: "#4b5563" }}>件まで</span>
            </div>
            <span className="seq-hint">既存の BLAST 設定（ローカルDB）を使用</span>
          </label>
          <div className="seq-label" style={{ flex: 1 }}>
            BLAST 実行先 / ローカルDB
            <div className="blast-backend-row" style={{ flexDirection: "row", gap: "1.2rem" }}>
              <label>
                <input
                  type="checkbox"
                  checked={blastUseLocal}
                  onChange={(e) => setBlastUseLocal(e.target.checked)}
                />{" "}
                ローカル BLAST+
              </label>
            </div>
            {blastUseLocal && (
              <>
                <div className="blast-backend-row" style={{ flexDirection: "row", gap: "1rem", flexWrap: "wrap" }}>
                  {localDbOptions.map((opt) => (
                    <label key={opt.value}>
                      <input
                        type="checkbox"
                        checked={selectedLocalDbs.includes(opt.value)}
                        onChange={() => toggleLocalDb(opt.value)}
                      />{" "}
                      {opt.label}
                    </label>
                  ))}
                </div>
                <div className="primer-row">
                  <input
                    className="seq-input"
                    type="text"
                    value={customLocalDb}
                    onChange={(e) => setCustomLocalDb(e.target.value)}
                    placeholder="追加の makeblastdb プレフィックス (任意)"
                  />
                </div>
                <div className="tag-row">
                  <span className="tag-label">選択中</span>
                  <code className="tag-db">
                    {effectiveLocalDbs.length ? effectiveLocalDbs.join(", ") : "-"}
                  </code>
                </div>
                <span className="seq-hint">
                  BLAST DB base: {DEFAULT_BLAST_DB_BASE} ／ num_threads:{" "}
                  {blastNumThreads != null ? blastNumThreads : "自動 (CPU に応じて最大24、複数DBは自動で割り当て)"}
                </span>
              </>
            )}
            <div className="primer-row" style={{ gap: "1rem", flexWrap: "wrap", marginTop: "0.5rem" }}>
              <label className="seq-label" style={{ maxWidth: "180px" }}>
                task
                <select
                  className="seq-input"
                  value={blastTask}
                  onChange={(e) => setBlastTask(e.target.value)}
                >
                  <option value="blastn-short">blastn-short</option>
                  <option value="blastn">blastn</option>
                  <option value="megablast">megablast</option>
                </select>
              </label>
              <label className="seq-label" style={{ maxWidth: "200px" }}>
                E-value
                <input
                  className="seq-input"
                  type="number"
                  step="any"
                  value={blastEvalue}
                  onChange={(e) => setBlastEvalue(Number(e.target.value) || 1e-5)}
                />
              </label>
              <label className="seq-label" style={{ maxWidth: "200px" }}>
                max_hsps
                <input
                  className="seq-input"
                  type="number"
                  min={1}
                  value={blastMaxHsps ?? ""}
                  onChange={(e) => {
                    const v = Number(e.target.value);
                    setBlastMaxHsps(Number.isNaN(v) ? null : Math.max(1, v));
                  }}
                  placeholder="BLASTデフォルト"
                />
              </label>
              <label className="seq-label" style={{ maxWidth: "220px" }}>
                num_threads
                <input
                  className="seq-input"
                  type="number"
                  min={1}
                  value={blastNumThreads ?? ""}
                  onChange={(e) => {
                    const v = Number(e.target.value);
                    setBlastNumThreads(Number.isNaN(v) ? null : Math.max(1, v));
                  }}
                  placeholder="自動"
                />
                <span className="seq-hint">未指定なら CPU に応じて自動</span>
              </label>
              <label className="seq-label" style={{ maxWidth: "220px" }}>
                local mode
                <span className="seq-hint">CPU（通常）</span>
              </label>
            </div>
          </div>
        </div>

        <div className="primer-row" style={{ justifyContent: "flex-end" }}>
          <button
            type="button"
            className="seq-button"
            onClick={() => setShowManualInputs((v) => !v)}
            style={{ background: showManualInputs ? "linear-gradient(135deg, #2563eb, #1d4ed8)" : undefined }}
          >
            {showManualInputs ? "手動入力を隠す" : "手動入力を表示"}
          </button>
        </div>

        {showManualInputs && (
          <>
            <label className="seq-label">
              配列をドラッグしてターゲットにする（開始/終了を自動設定）:
              <textarea
                className="seq-textarea"
                rows={4}
                value={normalizedSeq}
                onMouseUp={(e) => {
                  const start = e.currentTarget.selectionStart;
                  const end = e.currentTarget.selectionEnd;
                  if (start !== end) {
                    setSelectionRange({ start: start + 1, end });
                    setSubStart(start + 1);
                    setSubEnd(end);
                  }
                }}
                readOnly
              />
              {selectionRange ? (
                <p className="seq-hint">
                  選択範囲: {selectionRange.start}–{selectionRange.end} bp
                </p>
              ) : (
                <p className="seq-hint">テキストを選択すると、その範囲をターゲットに設定します。</p>
              )}
            </label>
            <label className="seq-label">
              エキソン範囲 (1-based, 例: 100-200,300-450):
              <input
                className="seq-input"
                type="text"
                value={exonRanges}
                onChange={(e) => setExonRanges(e.target.value)}
                placeholder="100-200,300-450"
              />
            </label>
            <label className="seq-label">
              CDS 範囲 (任意, 例: 150-900):
              <input
                className="seq-input"
                type="text"
                value={cdsRange}
                onChange={(e) => setCdsRange(e.target.value)}
                placeholder="150-900"
              />
            </label>

            <div className="primer-row">
              <label className="seq-label">
                サブ領域開始:
                <input
                  className="seq-input"
                  type="number"
                  value={subStart ?? ""}
                  onChange={(e) => {
                    const v = Number(e.target.value);
                    setSubStart(Number.isNaN(v) ? null : v);
                  }}
                  placeholder="例: 200"
                />
              </label>
              <label className="seq-label">
                サブ領域終了:
                <input
                  className="seq-input"
                  type="number"
                  value={subEnd ?? ""}
                  onChange={(e) => {
                    const v = Number(e.target.value);
                    setSubEnd(Number.isNaN(v) ? null : v);
                  }}
                  placeholder="例: 400"
                />
              </label>
            </div>

            <div className="primer-row">
              <label className="seq-label">
                product size range:
                <input
                  className="seq-input"
                  type="text"
                  value={productSizeRange}
                  onChange={(e) => setProductSizeRange(e.target.value)}
                />
              </label>
              <label className="seq-label">
                候補数:
                <input
                  className="seq-input"
                  type="number"
                  value={numReturn}
                  min={1}
                  max={20}
                  onChange={(e) => setNumReturn(Number(e.target.value) || 5)}
                />
              </label>
            </div>
          </>
        )}
        {infoMessage && <p className="seq-hint">{infoMessage}</p>}
        <button type="button" className="seq-button" onClick={handleDesign} disabled={designLoading}>
          {designLoading ? "設計中..." : "プライマーを設計する"}
        </button>
        {designError && <p className="seq-error">エラー: {designError}</p>}
      </div>

      <div className="exon-visual">
        <div className="exon-visual-grid">
          <div>
            <h4>CDS/エキソン ハイライト</h4>
            {normalizedSeq && (parsedExons.length > 0 || parsedCds.length > 0) ? (
              <div ref={fastaHighlightRef}>
                <FeatureSequenceView
                  sequence={normalizedSeq}
                  header={`>${genomeHeader || geneId || "sequence"}（CDS/エキソン を枠で表示）`}
                  exonRanges={featureViewExons}
                  cdsRanges={featureViewCds}
                  primerRanges={featureViewPrimerRanges}
                  highlightRange={guideRegion ? { start: guideRegion.start, end: guideRegion.end } : null}
                  blockLen={FEATURE_VIEW_BLOCK_LEN}
                  fontSize="0.92rem"
                />
              </div>
            ) : (
              <p className="seq-hint">ゲノム配列とエキソン/CDS 範囲をセットするとここにハイライトが表示されます。</p>
            )}
          </div>

          <div className="primer-sidebar">
            <h4>プライマー抜粋</h4>
            <div className="primer-tabs">
              <button
                type="button"
                className={`primer-tab-btn ${primerPatternTab === "pattern1" ? "is-active" : ""}`}
                onClick={() => setPrimerPatternTab("pattern1")}
              >
                パターン1
              </button>
              <button
                type="button"
                className={`primer-tab-btn ${primerPatternTab === "pattern2" ? "is-active" : ""}`}
                onClick={() => setPrimerPatternTab("pattern2")}
              >
                パターン2
              </button>
            </div>
            {currentPatternPrimers.length === 0 && (
              <p className="seq-hint">設計後にここにプライマーを表示します。</p>
            )}
            {currentPatternPrimers.length > 0 && (
              <div className="primer-sidebar-list">
                {currentPatternPrimers.map((p, idx) => (
                  <div
                    key={`${p.exonIndex}-${p.subIndex || 0}-${primerPatternTab}-${idx}`}
                    className="primer-sidebar-card"
                    onClick={() => {
                      if (!normalizedSeq) return;
                      const leftStart = p.leftStart ?? 1;
                      const leftEnd = p.leftLen ? leftStart + p.leftLen - 1 : leftStart;
                      const rightStart = p.rightStart ?? leftStart;
                      const rightEnd = p.rightLen ? rightStart + p.rightLen - 1 : rightStart;
                      let regionStart = Math.min(leftStart, rightStart);
                      let regionEnd = Math.max(leftEnd, rightEnd);
                      if (p.product && p.product > 0) {
                        regionEnd = regionStart + p.product - 1;
                      }
                      regionStart = Math.max(1, regionStart);
                      regionEnd = Math.min(normalizedSeq.length, regionEnd);
                      scrollToRegion(regionStart, regionEnd);
                    }}
                  >
                    <div className="primer-sidebar-title">
                      {p.exonEndIndex && p.exonEndIndex !== p.exonIndex
                        ? `Exon ${p.exonIndex}–${p.exonEndIndex}`
                        : `Exon ${p.exonIndex}`}
                      {p.subIndex ? ` sub ${p.subIndex}` : ""}
                      {p.product ? ` / ${p.product}bp` : ""}
                    </div>
                    <div className="primer-sidebar-seq">
                      <span className="primer-tag">F</span>
                      <code>{p.left}</code>
                    </div>
                    <div className="primer-sidebar-seq">
                      <span className="primer-tag">R</span>
                      <code>{p.right}</code>
                    </div>
                    {p.quality && (
                      <div className="seq-hint" style={{ marginTop: "0.2rem" }}>
                        品質ランク: {p.quality}
                      </div>
                    )}
                    {p.note && (
                      <div className="seq-hint" style={{ marginTop: "0.1rem", color: "#9d4700" }}>
                        {p.note}
                      </div>
                    )}
                    <div className="primer-marker-track">
                      {p.leftStart && p.leftLen && (
                        <div
                          className="primer-marker primer-marker-left"
                          style={{
                            left: `${((p.leftStart - 1) / normalizedSeq.length) * 100}%`,
                            width: `${(p.leftLen / normalizedSeq.length) * 100}%`,
                          }}
                          title={`F ${p.leftStart}-${p.leftStart + p.leftLen - 1}`}
                        />
                      )}
                      {p.rightStart && p.rightLen && (
                        <div
                          className="primer-marker primer-marker-right"
                          style={{
                            left: `${((p.rightStart - 1) / normalizedSeq.length) * 100}%`,
                            width: `${(p.rightLen / normalizedSeq.length) * 100}%`,
                          }}
                          title={`R ${p.rightStart}-${p.rightStart + p.rightLen - 1}`}
                        />
                      )}
                    </div>
                  </div>
                ))}
              </div>
            )}
          </div>
        </div>
        {designResult && (
          <>
            <div className="table-scroll">
              <table className="seq-table">
                <thead>
                  <tr>
                    <th>
                      <input
                        type="checkbox"
                        checked={
                          designResult
                            ? selectedPairIndices.length === designResult.candidates.length
                            : false
                        }
                        onChange={(e) => {
                          if (!designResult) return;
                          if (e.target.checked) {
                            setSelectedPairIndices(designResult.candidates.map((_, i) => i));
                            setSelectedPair(designResult.candidates.length ? 0 : null);
                          } else {
                            setSelectedPairIndices([]);
                            setSelectedPair(null);
                          }
                        }}
                        aria-label="select-all"
                      />
                    </th>
                    <th>左プライマー</th>
                    <th>右プライマー</th>
                    <th>産物長</th>
                    <th>ペナルティ</th>
                  </tr>
                </thead>
                <tbody>
                  {designResult.candidates.map((pair, idx) => {
                    const checked = selectedPairIndices.includes(idx);
                    return (
                      <tr
                        key={idx}
                        className={checked ? "primer-row-selected" : undefined}
                        onClick={() => {
                          if (checked) {
                            setSelectedPairIndices(selectedPairIndices.filter((i) => i !== idx));
                          } else {
                            setSelectedPairIndices([...selectedPairIndices, idx]);
                          }
                          setSelectedPair(idx);
                        }}
                      >
                        <td>
                          <input
                            type="checkbox"
                            checked={checked}
                            onChange={(e) => {
                              e.stopPropagation();
                              if (e.target.checked) {
                                setSelectedPairIndices([...selectedPairIndices, idx]);
                                setSelectedPair(idx);
                              } else {
                                setSelectedPairIndices(selectedPairIndices.filter((i) => i !== idx));
                                setSelectedPair(null);
                              }
                            }}
                          />
                        </td>
                        <td>{pair.left_sequence}</td>
                        <td>{pair.right_sequence}</td>
                        <td>{pair.product_size ?? "-"}</td>
                        <td>{pair.pair_penalty?.toFixed(2) ?? "-"}</td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
            {renderAmplicon()}
            <div className="primer-blast-block">
              <div className="primer-blast-header">
                <div>
                  <h3>Primer-BLAST（選択ペアの左右を並列 BLAST）</h3>
                  <p className="panel-hint">ローカルDBで実行</p>
                </div>
                <div className="primer-blast-actions">
                  <div className="seq-label">
                    <div className="blast-backend-row checklist-grid">
                      <span>ローカル DB（複数選択可）:</span>
                      {localDbOptions.map((opt) => (
                        <label key={opt.value}>
                          <input
                            type="checkbox"
                            checked={selectedLocalDbs.includes(opt.value)}
                            onChange={() => toggleLocalDb(opt.value)}
                          />{" "}
                          {opt.label}
                        </label>
                      ))}
                    </div>
                    <div className="primer-row">
                      <input
                        className="seq-input"
                        type="text"
                        value={customLocalDb}
                        onChange={(e) => setCustomLocalDb(e.target.value)}
                        placeholder="追加の makeblastdb prefix (任意)"
                      />
                    </div>
                    <div className="tag-row">
                      <span className="tag-label">選択中</span>
                      <code className="tag-db">
                        {effectiveLocalDbs.length ? effectiveLocalDbs.join(", ") : "-"}
                      </code>
                    </div>
                  </div>
                  <label className="seq-label">
                    実行先:
                    <div className="blast-backend-row checklist-grid">
                      <label>
                        <input
                          type="checkbox"
                          checked={blastUseLocal}
                          onChange={(e) => setBlastUseLocal(e.target.checked)}
                        />{" "}
                        ローカル BLAST+
                      </label>
                    </div>
                  </label>
                  <label className="seq-label">
                    最大ヒット数:
                    <input
                      className="seq-input"
                      type="number"
                      min={1}
                      max={50}
                      value={blastMaxHits}
                      onChange={(e) =>
                        setBlastMaxHits(Number.isNaN(Number(e.target.value)) ? 5 : Number(e.target.value))
                      }
                    />
                  </label>
                  <button type="button" className="seq-button" onClick={runPrimerBlast} disabled={blastLoading}>
                    {blastLoading ? "BLAST 実行中..." : "選択ペアを BLAST"}
                  </button>
                  {blastUseLocal && blastLoading ? (
                    <JobProgressCard
                      title="ローカル BLAST"
                      jobId={blastJobId}
                      job={blastJobInfo}
                      onCancel={blastJobId ? cancelBlastJob : null}
                      cancelDisabled={!blastJobId}
                    />
                  ) : null}
                  {blastError && <p className="seq-error">エラー: {blastError}</p>}
                </div>
              </div>
              <div className="primer-blast-results">
                {primerBlastResults.length > 0 && (
                  <div className="table-scroll" style={{ marginBottom: "0.5rem" }}>
                    <table className="seq-table">
                      <thead>
                        <tr>
                          <th>Pair</th>
                          <th>側</th>
                          <th>source</th>
                          <th>ヒット数</th>
                          <th>Top ID</th>
                          <th>%id</th>
                          <th>E-value</th>
                        </tr>
                      </thead>
                      <tbody>
                        {primerBlastResults.flatMap((res) => {
                          const rows: Array<{
                            key: string;
                            pair: number;
                            side: string;
                            name: string;
                            count: number;
                            topId: string;
                            topPident: string;
                            topE: string;
                          }> = [];
                          const summarize = (side: string, resp: BlastResponse | null) => {
                            if (!resp) return;
                            const grouped = resp.hits.reduce<Record<string, BlastHit[]>>((acc, h) => {
                              const src = h.source || side;
                              acc[src] = acc[src] || [];
                              acc[src].push(h);
                              return acc;
                            }, {});
                            Object.entries(grouped).forEach(([src, hits]) =>
                              rows.push({
                                key: `${res.pairIndex}-${side}-${src}`,
                                pair: res.pairIndex + 1,
                                side,
                                name: src,
                                count: hits.length,
                                topId: hits[0]?.sseqid.split(/\s+/)[0] ?? "-",
                                topPident: hits[0] ? hits[0].pident.toFixed(1) : "-",
                                topE: hits[0] ? hits[0].evalue.toExponential(2) : "-",
                              }),
                            );
                          };
                          summarize("Left", res.left);
                          summarize("Right", res.right);
                          return rows;
                        }).map((r) => (
                          <tr key={r.key}>
                            <td>{r.pair}</td>
                            <td>{r.side}</td>
                            <td>{r.name}</td>
                            <td>{r.count}</td>
                            <td>{r.topId}</td>
                            <td>{r.topPident}</td>
                            <td>{r.topE}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>
                )}
                {primerBlastResults.map((res) => (
                  <div key={res.pairIndex} className="primer-blast-pair">
                    <h4 className="panel-title">Pair {res.pairIndex + 1}</h4>
                    {renderBlastTable(res.left, "Left primer")}
                    {renderBlastTable(res.right, "Right primer")}
                  </div>
                ))}
                {primerBlastResults.length === 0 && !blastLoading && (
                  <p className="seq-hint">ペアを選んで「選択ペアを BLAST」を押すとここに結果が表示されます。</p>
                )}
              </div>
            </div>
          </>
        )}
        {exonPrimerPlans.length > 0 && (
          <div className="seq-result-block">
            <h3>自動プライマー（Exon/CDS）</h3>
            {batchStatus && (
              <div className="primer-row" style={{ alignItems: "center" }}>
                <p className="seq-hint" style={{ marginBottom: 0 }}>
                  自動設計の進行状況: {batchStatus.done}/{batchStatus.total}（
                  {batchStatus.total > 0 ? ((batchStatus.done / batchStatus.total) * 100).toFixed(0) : "0"}%）
                  {batchStatus.cancelled ? "（キャンセル済み）" : ""}
                </p>
                {!batchStatus.cancelled && batchStatus.done < batchStatus.total && (
                  <button
                    type="button"
                    className="seq-button secondary danger"
                    style={{ marginLeft: "0.75rem" }}
                    onClick={() => {
                      batchCancelRef.current = true;
                    }}
                  >
                    自動設計をキャンセル
                  </button>
                )}
              </div>
            )}
            {coverageStats && (
              <p className="seq-hint">
                現在のパターン ({primerPatternTab === "pattern1" ? "パターン1" : "パターン2"}) のカバレッジ:
                {coverageStats.cds && coverageStats.cds.total > 0 && (
                  <> CDS {coverageStats.cds.hit}/{coverageStats.cds.total} bp</>
                )}
                {coverageStats.exon && coverageStats.exon.total > 0 && (
                  <> ／ Exon {coverageStats.exon.hit}/{coverageStats.exon.total} bp</>
                )}
              </p>
            )}
            <div className="primer-row" style={{ marginBottom: "0.4rem" }}>
              <button
                type="button"
                className="seq-button secondary"
                onClick={() => {
                  const md = buildMarkdownReport();
                  if (!md) return;
                  const base =
                    geneId?.trim() ||
                    (mode === "sequence" ? "sequence_primers" : "exon_primers");
                  downloadMarkdown(md, base);
                }}
              >
                結果を Markdown として保存
              </button>
              <button
                type="button"
                className="seq-button secondary"
                onClick={() => {
                  const md = buildMarkdownReport();
                  if (!md) return;
                  openPrintViewForMarkdown(
                    md,
                    geneId
                      ? `Exon/CDS レポート (${geneId})`
                      : "Exon/CDS プライマー設計レポート",
                  );
                }}
              >
                印刷用ビューを開く（PDF 保存に利用）
              </button>
            </div>
            {currentPatternPrimers.length > 0 && (
              <div style={{ marginBottom: "0.6rem" }}>
                <p className="seq-hint">
                  Primer 逆引きタブの「プライマー配列」欄に貼り付けられる形式（2 行ずつで 1 ペア）:
                </p>
                <textarea
                  className="seq-textarea"
                  rows={Math.min(2 * currentPatternPrimers.length + 2, 12)}
                  readOnly
                  value={currentPatternPrimers
                    .map((p) => `${p.left}\n${p.right}`)
                    .join("\n")}
                />
              </div>
            )}
            <div className="table-scroll" style={{ marginBottom: "0.6rem" }}>
              <table className="seq-table">
                <thead>
                  <tr>
                    <th>Exon</th>
                    <th>範囲</th>
                    <th>サブ</th>
                    <th>候補数</th>
                    <th>Status</th>
                  </tr>
                </thead>
                <tbody>
                  {exonPrimerPlans.map((p) => (
                    <tr key={`${p.exonIndex}-${p.subIndex || 0}`}>
                      <td>{p.exonIndex}</td>
                      <td>
                        {p.range[0]}–{p.range[1]}
                      </td>
                      <td>{p.subIndex ?? (p.exonEndIndex && p.exonEndIndex !== p.exonIndex
                        ? `Exon${p.exonIndex}–${p.exonEndIndex}`
                        : "-")}</td>
                      <td>{p.result?.candidates.length ?? 0}</td>
                      <td>
                        {p.loading
                          ? "設計中"
                          : p.error
                            ? "エラー"
                            : p.result
                              ? "完了"
                              : "-"}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            {topCandidateSummary.length > 0 && (
              <div className="table-scroll" style={{ marginBottom: "0.6rem" }}>
                <table className="seq-table">
                  <thead>
                    <tr>
                      <th>ターゲット</th>
                      <th>範囲</th>
                      <th>左プライマー</th>
                      <th>右プライマー</th>
                      <th>産物長</th>
                      <th>設計スコア (S〜D, penalty)</th>
                    </tr>
                  </thead>
                  <tbody>
                    {topCandidateSummary.map((c) => (
                      <tr key={`${c.exonIndex}-${c.range[0]}-${c.range[1]}`}>
                        <td>{c.kind === "cds" ? `CDS${c.exonIndex}` : `Exon${c.exonIndex}`}</td>
                        <td>
                          {c.range[0]}–{c.range[1]}
                        </td>
                        <td>{c.left}</td>
                        <td>{c.right}</td>
                        <td>{c.product ?? "-"}</td>
                        <td>
                          {c.quality && (
                            <span
                              className={`quality-badge grade-${c.quality.toLowerCase()}`}
                              style={{ marginRight: "0.25rem" }}
                            >
                              {c.quality}
                            </span>
                          )}
                          {typeof c.penalty === "number"
                            ? c.penalty.toFixed(2)
                            : c.quality
                              ? ""
                              : "-"}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            )}
            {exonPrimerPlans.map((plan) => (
              <div key={`${plan.exonIndex}-${plan.subIndex || 0}`} className="primer-exon-block">
                <h4>
                  {(plan.kind ?? "exon") === "cds"
                    ? `CDS ${plan.exonIndex}`
                    : `Exon ${plan.exonIndex}`}{" "}
                  ({plan.range[0]}–{plan.range[1]})
                </h4>
                {plan.loading && <p className="seq-hint">設計中...</p>}
                {plan.error && <p className="seq-error">エラー: {plan.error}</p>}
                {plan.result && (
                  <div className="table-scroll">
                    <table className="seq-table">
                      <thead>
                        <tr>
                          <th>#</th>
                          <th>左プライマー</th>
                          <th>右プライマー</th>
                          <th>産物長</th>
                          <th>設計スコア (S〜D, penalty)</th>
                          <th>位置 / BLAST / 品質</th>
                        </tr>
                      </thead>
                      <tbody>
                        {plan.result.candidates.length === 0 && (
                          <tr>
                            <td colSpan={6}>
                              <span className="seq-hint">候補が見つかりませんでした。</span>
                            </td>
                          </tr>
                        )}
                        {(plan.filteredCandidates || plan.result.candidates).map((pair, idx) => {
                          const p = (pair as any).pair ? (pair as any).pair : (pair as PrimerPair);
                          const left = (pair as any).left as BlastResponse | undefined;
                          const right = (pair as any).right as BlastResponse | undefined;
                          const note = (pair as any).note as string | undefined;
                          const ampCount =
                            (pair as any).ampliconCount as number | undefined;
                          const quality =
                            (pair as any).quality as "S" | "A" | "B" | "C" | "D" | undefined;
                          return (
                            <tr key={idx}>
                              <td>{idx + 1}</td>
                              <td>{p.left_sequence}</td>
                              <td>{p.right_sequence}</td>
                              <td>{p.product_size ?? "-"}</td>
                              <td>
                                {quality && (
                                  <span
                                    className={`quality-badge grade-${quality.toLowerCase()}`}
                                    style={{ marginRight: "0.25rem" }}
                                  >
                                    {quality}
                                  </span>
                                )}
                                {typeof p.pair_penalty === "number"
                                  ? p.pair_penalty.toFixed(2)
                                  : quality
                                    ? ""
                                    : "-"}
                              </td>
                              <td>
                                L:{p.left_start} / R:{p.right_start}
                                {left && right && (
                                  <>
                                    <div className="seq-hint">
                                      BLAST ヒット数 (合計): L {left.num_hits ?? 0} / R{" "}
                                      {right.num_hits ?? 0}
                                    </div>
                                    {(() => {
                                      const perDb = summarizeBlastHitsPerDb(left, right);
                                      return perDb.length ? (
                                        <div className="seq-hint">
                                          DB別ヒット数:{" "}
                                          {perDb.map((d, idx) => (
                                            <span key={d.db}>
                                              {idx > 0 && " / "}
                                              {d.db}: L {d.left} / R {d.right}
                                            </span>
                                          ))}
                                        </div>
                                      ) : null;
                                    })()}
                                  </>
                                )}
                                {typeof ampCount === "number" && (
                                  <div className="seq-hint">
                                    予測PCR産物（ローカル）: {ampCount} 本
                                  </div>
                                )}
                                {note && (
                                  <div className="seq-hint" style={{ color: "#9d4700" }}>
                                    {note}
                                  </div>
                                )}
                                <div className="primer-row" style={{ marginTop: "0.3rem", gap: "0.25rem", flexWrap: "wrap" }}>
                                  {setPresetReversePair && setActiveTab && (
                                    <button
                                      type="button"
                                      className="seq-button secondary"
                                      style={{ padding: "0.15rem 0.4rem", fontSize: "0.75rem" }}
                                      onClick={() => {
                                        setPresetReversePair({
                                          primer1: p.left_sequence,
                                          primer2: p.right_sequence,
                                        });
                                        setActiveTab("primer_reverse");
                                      }}
                                    >
                                      逆引きタブで評価
                                    </button>
                                  )}
                                  <button
                                    type="button"
                                    className="seq-button secondary"
                                    style={{ padding: "0.15rem 0.4rem", fontSize: "0.75rem" }}
                                    onClick={() => {
                                      const base = new URL(window.location.href);
                                      base.pathname = "/";
                                      base.hash = "";
                                      base.port = "3100";
                                      const qp = new URLSearchParams();
                                      qp.set("f", p.left_sequence);
                                      qp.set("r", p.right_sequence);
                                      base.search = qp.toString();
                                      window.open(base.toString(), "_blank", "noopener,noreferrer");
                                    }}
                                  >
                                    AB1 Viewerで開く
                                  </button>
                                </div>
                              </td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                    <p className="seq-hint">
                      {mode === "sequence" || sequencingMode
                        ? `シーケンスプライマーとして、産物長を ${seqProductMin}–${seqProductMax} bp に設定して設計しています。`
                        : `ターゲット exon 長 ${plan.range[1] - plan.range[0] + 1} bp に合わせ、産物長を広めに設定しています (約${Math.max(80, plan.range[1] - plan.range[0] + 1)}–${Math.max(80, plan.range[1] - plan.range[0] + 1) + 150} bp)。`}
                    </p>
                  </div>
                )}
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
};

function parseRanges(text: string): Array<[number, number]> {
  return text
    .split(",")
    .map((s) => s.trim())
    .filter(Boolean)
    .map((part) => {
      const [a, b] = part.split("-").map((x) => Number(x));
      if (!Number.isFinite(a) || !Number.isFinite(b)) return null;
      const start = Math.min(a, b);
      const end = Math.max(a, b);
      return [start, end] as [number, number];
    })
    .filter((r): r is [number, number] => !!r);
}

// Synthetic sample data (not derived from any organism).
const sampleCds = "ATG" + "GCT".repeat(200) + "TAA";
const sampleExons = "ATG" + "GCT".repeat(100) + "\n" + "GCT".repeat(100) + "TAA";





