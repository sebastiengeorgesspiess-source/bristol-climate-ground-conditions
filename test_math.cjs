const assert=require('assert'),M=require('./static/bristol-math-v51.js');
assert.deepEqual(M.rank([2,2,8]),[1.5,1.5,3]);assert.equal(M.corr([1,1,1],[1,2,3]),null);
assert.deepEqual(M.transform([{date:'2025-01',value:120}],'z'),[]);
assert.deepEqual(M.smooth([{date:'2025-01',value:1},{date:'2025-03',value:3}],2),[]);
const q=M.bh([{p:.001},{p:.4},{p:.8}]);assert(q[0]<.05&&q[1]>.05);
assert.equal(M.infer([{date:'2025-01',a:1,b:2}]),null);
const a=[{date:'2025-12',value:1}],b=[{date:'2026-01',value:2}];assert.equal(M.pairs(a,b,'raw',1,1,2025,2025).length,0);
console.log('PASS: tied ranks, constant/short series, missing baseline, calendar gaps, selected BH q, date-boundary lag');
