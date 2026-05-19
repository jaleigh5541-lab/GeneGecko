import React from "react";

interface Props {
  onClick: () => void;
  isProcessing: boolean;
  loadingText: string;
}

export default function ProcessButton({ onClick, isProcessing, loadingText }: Props) {
  return (
    <button className="process-btn" onClick={onClick} disabled={isProcessing}>
      {isProcessing ? (
        <>
          <span className="gecko-steps">
            <span className="gecko-toe" />
            <span className="gecko-toe" />
            <span className="gecko-toe" />
            <span className="gecko-toe" />
            <span className="gecko-toe" />
          </span>
          {loadingText}
        </>
      ) : (
        "Go Gecko!"
      )}
    </button>
  );
}
