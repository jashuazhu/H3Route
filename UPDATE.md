# Updating Content

This project converts `webpage.docx` into a styled `index.html` using a minimal PowerShell-based DOCX reader (no external tools).

Regenerate after editing `webpage.docx`:

```powershell
Add-Type -AssemblyName System.IO.Compression, System.IO.Compression.FileSystem
$ErrorActionPreference = "Stop"
# Use the conversion block used to generate index.html (embedded in automation).
```

Notes:
- The converter preserves paragraphs, simple bold/italic, basic ordered lists, and the first two centered lines as title/subtitle.
- Images, tables, complex numbering styles, and hyperlinks may not be fully preserved.
- `CONTENT.md` is a rough Markdown export intended for quick diffing or reuse.

Styling: `styles.css` defines the book-like layout applied to `index.html`.
