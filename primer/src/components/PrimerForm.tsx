import { useState } from 'react';

interface PrimerOptions {
  addFwd: string;
  addRev: string;
  addFwdLen: [number, number];
  addRevLen: [number, number];
  offtargetCheck: string;
  optimalTm: number;
  optimalGc: number;
  optimalLen: number;
  penaltyTm: number;
  penaltyGc: number;
  penaltyLen: number;
  penaltyTmDiff: number;
  penaltyDg: number;
  penaltyOffTarget: number;
}

const defaultOptions: PrimerOptions = {
  addFwd: '',
  addRev: '',
  addFwdLen: [-1, -1],
  addRevLen: [-1, -1],
  offtargetCheck: '',
  optimalTm: 62.0,
  optimalGc: 0.5,
  optimalLen: 22,
  penaltyTm: 1.0,
  penaltyGc: 0.2,
  penaltyLen: 0.5,
  penaltyTmDiff: 1.0,
  penaltyDg: 2.0,
  penaltyOffTarget: 20.0,
};

// eslint-disable-next-line @typescript-eslint/no-explicit-any
type SubmitHandler = (data: any) => void;

interface PrimerFormProps {
  onSubmit: SubmitHandler;
  mode?: 'create' | 'score';
}

export default function PrimerForm({ onSubmit, mode = 'create' }: PrimerFormProps) {
  const [sequence, setSequence] = useState<string>('');
  const [fwdPrimer, setFwdPrimer] = useState<string>('');
  const [revPrimer, setRevPrimer] = useState<string>('');
  const [showAdvanced, setShowAdvanced] = useState<boolean>(false);
  const [options, setOptions] = useState<PrimerOptions>(defaultOptions);

  const handleSubmit = (e: React.FormEvent<HTMLFormElement>): void => {
    e.preventDefault();
    if (mode === 'create') {
      onSubmit({ sequence, options });
    } else {
      onSubmit({ fwdPrimer, revPrimer, sequence, options });
    }
  };

  const handleOptionChange = (key: keyof PrimerOptions, value: string | number): void => {
    setOptions((prev) => ({ ...prev, [key]: value }));
  };

  const sanitizeDNA = (value: string): string => {
    return value.toUpperCase().replace(/[^ATGC]/gi, '');
  };

  return (
    <form onSubmit={handleSubmit} className="primer-form">
      {mode === 'create' ? (
        <div className="form-group">
          <label htmlFor="sequence">DNA Sequence to Amplify</label>
          <textarea
            id="sequence"
            value={sequence}
            onChange={(e) => setSequence(sanitizeDNA(e.target.value))}
            placeholder="Enter DNA sequence (ATGC only)..."
            rows={5}
            required
          />
          <small>
            Enter the template sequence you want to amplify with PCR primers.
          </small>
        </div>
      ) : (
        <>
          <div className="form-group">
            <label htmlFor="fwdPrimer">Forward Primer</label>
            <input
              type="text"
              id="fwdPrimer"
              value={fwdPrimer}
              onChange={(e) => setFwdPrimer(sanitizeDNA(e.target.value))}
              placeholder="Enter forward primer sequence..."
              required
            />
          </div>
          <div className="form-group">
            <label htmlFor="revPrimer">Reverse Primer (optional)</label>
            <input
              type="text"
              id="revPrimer"
              value={revPrimer}
              onChange={(e) => setRevPrimer(sanitizeDNA(e.target.value))}
              placeholder="Enter reverse primer sequence..."
            />
          </div>
          <div className="form-group">
            <label htmlFor="seqScore">Template Sequence (optional)</label>
            <textarea
              id="seqScore"
              value={sequence}
              onChange={(e) => setSequence(sanitizeDNA(e.target.value))}
              placeholder="Enter template sequence for off-target checking..."
              rows={3}
            />
          </div>
        </>
      )}

      <div className="form-group">
        <button
          type="button"
          className="toggle-advanced"
          onClick={() => setShowAdvanced(!showAdvanced)}
        >
          {showAdvanced ? 'Hide' : 'Show'} Advanced Options
        </button>
      </div>

      {showAdvanced && (
        <div className="advanced-options">
          {mode === 'create' && (
            <>
              <div className="form-row">
                <div className="form-group">
                  <label htmlFor="addFwd">Add to FWD Primer (5&apos; to 3&apos;)</label>
                  <input
                    type="text"
                    id="addFwd"
                    value={options.addFwd}
                    onChange={(e) => handleOptionChange('addFwd', sanitizeDNA(e.target.value))}
                    placeholder="e.g., GGTCTC for BsaI"
                  />
                </div>
                <div className="form-group">
                  <label htmlFor="addRev">Add to REV Primer (5&apos; to 3&apos;)</label>
                  <input
                    type="text"
                    id="addRev"
                    value={options.addRev}
                    onChange={(e) => handleOptionChange('addRev', sanitizeDNA(e.target.value))}
                    placeholder="e.g., GAAGAC for BpiI"
                  />
                </div>
              </div>
            </>
          )}

          <div className="form-group">
            <label htmlFor="offtargetCheck">Off-target Check Sequence</label>
            <textarea
              id="offtargetCheck"
              value={options.offtargetCheck}
              onChange={(e) => handleOptionChange('offtargetCheck', sanitizeDNA(e.target.value))}
              placeholder="Optional: Enter sequence to check for off-target binding sites..."
              rows={2}
            />
          </div>

          <h4>Optimal Parameters</h4>
          <div className="form-row">
            <div className="form-group">
              <label htmlFor="optimalTm">Optimal Tm (&deg;C)</label>
              <input
                type="number"
                id="optimalTm"
                value={options.optimalTm}
                onChange={(e) => handleOptionChange('optimalTm', parseFloat(e.target.value))}
                step="0.1"
              />
            </div>
            <div className="form-group">
              <label htmlFor="optimalGc">Optimal GC Ratio</label>
              <input
                type="number"
                id="optimalGc"
                value={options.optimalGc}
                onChange={(e) => handleOptionChange('optimalGc', parseFloat(e.target.value))}
                step="0.01"
                min="0"
                max="1"
              />
            </div>
            <div className="form-group">
              <label htmlFor="optimalLen">Optimal Length (bp)</label>
              <input
                type="number"
                id="optimalLen"
                value={options.optimalLen}
                onChange={(e) => handleOptionChange('optimalLen', parseInt(e.target.value, 10))}
              />
            </div>
          </div>

          <h4>Penalty Weights</h4>
          <div className="form-row">
            <div className="form-group">
              <label htmlFor="penaltyTm">Tm Penalty</label>
              <input
                type="number"
                id="penaltyTm"
                value={options.penaltyTm}
                onChange={(e) => handleOptionChange('penaltyTm', parseFloat(e.target.value))}
                step="0.1"
              />
            </div>
            <div className="form-group">
              <label htmlFor="penaltyGc">GC Penalty</label>
              <input
                type="number"
                id="penaltyGc"
                value={options.penaltyGc}
                onChange={(e) => handleOptionChange('penaltyGc', parseFloat(e.target.value))}
                step="0.1"
              />
            </div>
            <div className="form-group">
              <label htmlFor="penaltyLen">Length Penalty</label>
              <input
                type="number"
                id="penaltyLen"
                value={options.penaltyLen}
                onChange={(e) => handleOptionChange('penaltyLen', parseFloat(e.target.value))}
                step="0.1"
              />
            </div>
          </div>
          <div className="form-row">
            <div className="form-group">
              <label htmlFor="penaltyTmDiff">Tm Diff Penalty</label>
              <input
                type="number"
                id="penaltyTmDiff"
                value={options.penaltyTmDiff}
                onChange={(e) => handleOptionChange('penaltyTmDiff', parseFloat(e.target.value))}
                step="0.1"
              />
            </div>
            <div className="form-group">
              <label htmlFor="penaltyDg">dG Penalty</label>
              <input
                type="number"
                id="penaltyDg"
                value={options.penaltyDg}
                onChange={(e) => handleOptionChange('penaltyDg', parseFloat(e.target.value))}
                step="0.1"
              />
            </div>
            <div className="form-group">
              <label htmlFor="penaltyOffTarget">Off-target Penalty</label>
              <input
                type="number"
                id="penaltyOffTarget"
                value={options.penaltyOffTarget}
                onChange={(e) => handleOptionChange('penaltyOffTarget', parseFloat(e.target.value))}
                step="1"
              />
            </div>
          </div>
        </div>
      )}

      <button type="submit" className="submit-btn">
        {mode === 'create' ? 'Design Primers' : 'Score Primers'}
      </button>
    </form>
  );
}
