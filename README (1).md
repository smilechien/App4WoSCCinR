# WoS Stable Pipeline (Routine 1/2/3 + Summary)

## Run


 Chien, T. W., Yang, S., Chou, W., & Wang, W. C. (2026). App4WoSCCinR: A Web-Based Bibliometric Analysis Platform for the Web of Science Core Collection. https://smilechien.shinyapps.io/woscc/
 
🔐 CMC (Pass Code) Requirements

The CMC is a 10-digit pass code, available upon request from the author.

The "test" CMC provides a 7-day trial access.

🔐 In Chinese
CMC是10碼Pass Code(需向作者取用)
CMC的一次上載Data試用
Demo 不受限制
A.上載WoSCC下載來的.xls檣為限
   再按Run (WoS.....)
B. 網址是帶參數進來的URL, 仍需CMC
    (通常是帶入URL, 當****CMC****生效時, 直接呈現其結果,往下拉,百張視覺圖)
C. Demo 也需CMC

📥 Data Loading Instructions
A. Upload Local File

Upload a WoS Core Collection (.xls) file downloaded from WoSCC.

After uploading, click “Run (WoS → wide32/long32/long16 → metatable)” to process.

Only .xls files exported from WoSCC are supported for upload.

B. Load via URL

The URL option supports loading a file with parameters passed directly in the link.

A valid CMC is still required, even when using URL parameters.

C. Demo Mode

The Demo function also requires a valid CMC.

Trial access is available using test (7 days).



## Outputs
- `woswide32.csv`
- `woslong32.csv`
- `woslong16.csv`
- `wosmetatable.csv` (A1..A13 from Keywords Plus)
- `summary_report.csv`
- `summary_report.png`
