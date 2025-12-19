import { FC } from 'react';

interface PiecewiseScores {
  tm?: number;
  gc?: number;
  terminal3DG?: number;
  offTarget?: number;
  hairpin?: number;
  homodimer?: number;
  heterodimer?: number;
  tmDiff?: number;
  [key: string]: number | undefined;
}

interface QualityInfo {
  tier?: string;
  label?: string;
}

interface ScoreBreakdownPopupProps {
  compositeScore: number;
  rawScore?: number;
  criticalWarnings?: number;
  quality: QualityInfo;
  forwardScores?: PiecewiseScores;
  reverseScores?: PiecewiseScores;
  hasTemplate?: boolean;
  onClose: () => void;
}

declare const ScoreBreakdownPopup: FC<ScoreBreakdownPopupProps>;
export default ScoreBreakdownPopup;
