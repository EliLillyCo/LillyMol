// Exercise the embedded script with a small DOM harness; no web dependencies.
// Usage: node test/medchem_transform_activity_html_ui_test.js /path/to/report.html
'use strict';
const fs=require('fs'),vm=require('vm'),assert=require('assert');
const page=fs.readFileSync(process.argv[2],'utf8');
const data=JSON.parse(page.match(/<script id="data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
const script=page.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
class Element{
 constructor(tag){this.tag=tag;this.children=[];this.dataset={};this.attrs={};this.events={};this.className='';this.value='';this.classList={add:name=>{this.className+=' '+name;}};}
 append(...nodes){for(const node of nodes){if(node.tag==='fragment')this.children.push(...node.children);else this.children.push(node);}}
 replaceChildren(...nodes){this.children=[];this.append(...nodes);}
 setAttribute(k,v){this.attrs[k]=v;}
 addEventListener(k,f){this.events[k]=f;}
 click(){this.clicked=true;}
}
const elements={};for(const id of ['search','minObs','minPred','match','source','metric','direction','export','explanation','rows','empty','status','data'])elements[id]=new Element('div');
Object.assign(elements.data,{textContent:JSON.stringify(data)});
for(const [id,value] of Object.entries({search:'',minObs:'0',minPred:'0',match:'all',source:'observed',metric:'median',direction:'-1'}))elements[id].value=value;
let capturedBlob;
const document={getElementById:id=>elements[id],createElement:tag=>new Element(tag),createElementNS:(_ns,tag)=>new Element(tag),createDocumentFragment:()=>new Element('fragment')};
const context=vm.createContext({document,Blob,URL:{createObjectURL:blob=>{capturedBlob=blob;return 'blob:test';},revokeObjectURL:()=>{}},setTimeout:f=>f(),console});
vm.runInContext(script,context);
const names=()=>Array.from(vm.runInContext('visible.map(r=>r.name)',context));
function change(id,value){elements[id].value=String(value);elements[id].events.input();}
const expected=data.rows.slice().sort((a,b)=>(b.observed.median??-Infinity)-(a.observed.median??-Infinity)||(a.name<b.name?-1:1)).map(r=>r.name);
assert.deepEqual(names(),expected);
change('direction','1');const known=expected.filter(name=>data.rows.find(r=>r.name===name).observed.n),missing=expected.filter(name=>!data.rows.find(r=>r.name===name).observed.n);assert.deepEqual(names(),known.reverse().concat(missing));
change('minObs',999999);assert.equal(names().length,0);assert.equal(elements.empty.hidden,false);
change('minPred',1);change('match','any');assert.equal(names().length,data.rows.filter(r=>r.predicted.n>=1).length);
change('minObs',0);change('minPred',0);change('search',data.rows[0].name);assert(names().includes(data.rows[0].name));
change('search','');change('source','predicted');change('metric','mean');change('direction','-1');
const values=Array.from(vm.runInContext('visible.filter(r=>r.predicted.n).map(r=>r.predicted.mean)',context));
assert(values.every((v,i)=>i===0||values[i-1]>=v));
const checkbox=elements.rows.children[0].children[0].children[0].children[0];checkbox.checked=true;checkbox.onchange();
assert(elements.status.textContent.includes('1 shortlisted'));
change('search','this matches nothing');assert(elements.status.textContent.includes('1 shortlisted'));
change('search','');
let details;
function walk(e){if(e.tag==='details')details=e;for(const c of e.children||[])walk(c);}
walk(elements.rows);
if(details){const before=details.children.length;details.open=true;details.events.toggle();assert(details.children.length>before);assert(details.children.length<=data.max_examples);const after=details.children.length;details.events.toggle();assert.equal(details.children.length,after);}
elements.export.onclick();
capturedBlob.text().then(text=>{assert(text.startsWith('"transformation",'));assert(text.includes('\r\n'));assert.equal(text.trim().split('\r\n').length,2);console.log('HTML sorting, filtering, capped expansion, shortlist persistence and CSV export passed');}).catch(e=>{console.error(e);process.exitCode=1;});
