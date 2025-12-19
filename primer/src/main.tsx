import React from 'react';
import ReactDOM from 'react-dom/client';
import App from './App';
import './styles/index.css';      // Tailwind + custom components
import './styles/seqviz.css';     // SeqViz library styles (do not convert)
import './App.css';               // Legacy CSS (gradually migrate to Tailwind)

const rootElement = document.getElementById('root');

if (!rootElement) {
  throw new Error('Failed to find the root element');
}

ReactDOM.createRoot(rootElement).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
