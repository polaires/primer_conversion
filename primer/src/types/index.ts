/**
 * Core Type Definitions for Primer Design Tool
 *
 * These types cover the main data structures used throughout the application
 * for primer design, scoring, assembly, and visualization.
 */

// ============================================================================
// Primer Types
// ============================================================================

export interface Primer {
  sequence: string;
  tm: number;
  gc: number;
  length: number;
  direction?: 'forward' | 'reverse';
  score?: number;
  penalty?: number;
  hairpin?: HairpinResult;
  homodimer?: DimerResult;
  position?: number;
  bind?: [number, number];
  terminal3dG?: number;
  warnings?: string[];
}

export interface PrimerPair {
  forward: Primer;
  reverse: Primer;
  tmDifference: number;
  heterodimer?: DimerResult;
  quality: QualityTier;
}

export type QualityTier = 'excellent' | 'good' | 'acceptable' | 'poor';

// ============================================================================
// Thermodynamics Types
// ============================================================================

export interface TmResult {
  tm: number;
  method: 'Q5' | 'SantaLucia' | 'DNA24';
  parameters?: TmParameters;
}

export interface TmParameters {
  primerConc: number;    // nM
  saltConc: number;      // mM
  mgConc: number;        // mM
  dNTPConc: number;      // mM
}

export interface HairpinResult {
  dG: number;
  tm: number;
  structure?: string;
  position?: number;
  score?: number;
}

export interface DimerResult {
  dG: number;
  tm: number;
  position?: number;
  severity: 'none' | 'low' | 'moderate' | 'high';
  score?: number;
}

// ============================================================================
// Scoring Types
// ============================================================================

export interface ScoringWeights {
  tm?: number;
  gc?: number;
  terminal3dG?: number;
  tmDiff?: number;
  hairpin?: number;
  homodimer?: number;
  heterodimer?: number;
  offTarget?: number;
  length?: number;
  gcClamp?: number;
  homopolymer?: number;
  threeprimeComp?: number;
}

export interface ScoreBreakdown {
  total: number;
  components: {
    [key: string]: {
      score: number;
      weight: number;
      contribution: number;
      details?: string;
    };
  };
}

export type ScoringPreset = 'standard' | 'mutagenesis' | 'sequencing' | 'assembly' | 'strict';

// ============================================================================
// Assembly Types
// ============================================================================

export type AssemblyMethod = 'golden_gate' | 'gibson' | 'nebuilder_hifi';

export interface AssemblyFragment {
  id: string;
  name: string;
  sequence: string;
  length: number;
  type?: FragmentType;
  annotations?: Annotation[];
}

export type FragmentType = 'promoter' | 'rbs' | 'cds' | 'terminator' | 'vector' | 'insert' | 'other';

export interface GoldenGateEnzyme {
  name: string;
  recognition: string;
  overhang: number;
  temperature: number;
  cutSite?: number;
  methylationSensitive?: boolean;
}

export interface AssemblyResult {
  success: boolean;
  fragments: AssemblyFragment[];
  primers?: PrimerPair[];
  overhangs?: string[];
  fidelity?: number;
  warnings?: string[];
}

export interface Overhang {
  sequence: string;
  position: number;
  fidelity?: number;
  crossReactivity?: number;
}

// ============================================================================
// Mutagenesis Types
// ============================================================================

export type MutationType = 'substitution' | 'insertion' | 'deletion' | 'codon_change';

export interface MutationSpec {
  type: MutationType;
  position: number;
  original?: string;
  replacement?: string;
  codon?: string;
  aminoAcid?: string;
}

export interface MutagenesisResult {
  forward: Primer;
  reverse: Primer;
  mutantSequence: string;
  mutationType: MutationType;
  verified: boolean;
  warnings?: string[];
}

// ============================================================================
// Sequencing Types
// ============================================================================

export interface SequencingPrimer extends Primer {
  position: number;
  coverage: [number, number];  // Start and end positions covered
  direction: 'forward' | 'reverse';
}

export interface SequencingDesign {
  primers: SequencingPrimer[];
  coverage: number;  // Percentage
  gaps?: [number, number][];
  totalReads: number;
}

// ============================================================================
// Annotation Types
// ============================================================================

export interface Annotation {
  id?: string;
  name: string;
  type: AnnotationType;
  start: number;
  end: number;
  direction?: 1 | -1;
  color?: string;
  notes?: string;
}

export type AnnotationType =
  | 'gene'
  | 'CDS'
  | 'promoter'
  | 'terminator'
  | 'RBS'
  | 'primer_bind'
  | 'misc_feature'
  | 'restriction_site'
  | 'origin'
  | 'rep_origin';

// ============================================================================
// UI State Types
// ============================================================================

export type ToolMode =
  | 'primer-designer'
  | 'sequencing'
  | 'score'
  | 'tm'
  | 'viewer'
  | 'assembly-studio'
  | 'fragment-planner';

export type ToolCategory = 'design' | 'analysis' | 'assembly';

export interface Tool {
  id: ToolMode;
  label: string;
  category: ToolCategory;
  description: string;
}

// ============================================================================
// Design Options Types
// ============================================================================

export interface PrimerDesignOptions {
  minTm?: number;
  maxTm?: number;
  optimalTm?: number;
  minGC?: number;
  maxGC?: number;
  minLength?: number;
  maxLength?: number;
  circular?: boolean;
  offtargetCheck?: boolean;
  tmParams?: TmParameters;
}

export interface SequencingDesignOptions {
  readLength?: number;
  minCoverage?: number;
  overlap?: number;
  direction?: 'forward' | 'reverse' | 'both';
}

export interface MutagenesisDesignOptions {
  flankingLength?: number;
  minTm?: number;
  maxTm?: number;
  maxPrimerLength?: number;
  method?: 'quikchange' | 'overlap' | 'neb_q5';
}

// ============================================================================
// File Format Types
// ============================================================================

export interface GenBankFeature {
  type: string;
  location: string;
  qualifiers: Record<string, string[]>;
}

export interface GenBankRecord {
  name: string;
  sequence: string;
  features: GenBankFeature[];
  annotations?: Record<string, string>;
  circular?: boolean;
  length: number;
}

export interface FastaRecord {
  id: string;
  description?: string;
  sequence: string;
}

// ============================================================================
// Utility Types
// ============================================================================

export type Direction = 'forward' | 'reverse';

export interface Range {
  start: number;
  end: number;
}

export interface Position {
  index: number;
  base: string;
}

// React Event Handler Types (convenience exports)
export type InputChangeHandler = React.ChangeEvent<HTMLInputElement>;
export type TextAreaChangeHandler = React.ChangeEvent<HTMLTextAreaElement>;
export type SelectChangeHandler = React.ChangeEvent<HTMLSelectElement>;
export type FormSubmitHandler = React.FormEvent<HTMLFormElement>;
export type ButtonClickHandler = React.MouseEvent<HTMLButtonElement>;
