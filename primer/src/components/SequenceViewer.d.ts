import { FC } from 'react';

interface Annotation {
  name: string;
  start: number;
  end: number;
  color: string;
  direction: number;
  type?: string;
}

interface HighlightRegion {
  start: number;
  end: number;
  color?: string;
}

interface SequenceViewerProps {
  sequence: string;
  name?: string;
  fwdPrimer?: string;
  revPrimer?: string;
  addFwd?: string;
  addRev?: string;
  viewer?: 'linear' | 'circular' | 'both';
  height?: number;
  annotations?: Annotation[];
  showControls?: boolean;
  showInfo?: boolean;
  showLegend?: boolean;
  compact?: boolean;
  enzymes?: string[];
  highlightRegions?: HighlightRegion[];
  circular?: boolean;
  onSequenceChange?: (sequence: string) => void;
  editable?: boolean;
}

interface AssemblyViewerProps {
  fragments: Array<{
    sequence: string;
    name: string;
    annotations?: Annotation[];
  }>;
  assemblyType?: 'golden-gate' | 'gibson' | 'restriction' | 'other';
  height?: number;
}

declare const SequenceViewer: FC<SequenceViewerProps>;
export const AssemblyViewer: FC<AssemblyViewerProps>;
export default SequenceViewer;
