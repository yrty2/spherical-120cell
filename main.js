const canvas=document.querySelector(".canvas");
let cam=[0,0,0];
let rot=[1,0,0,0];
let vert=[];
function rect(s){
    return [[s,s,s],[s,s,-s],[s,-s,s],[-s,s,s],[-s,-s,s],[-s,s,-s],[s,-s,-s],[-s,-s,-s]];
}    
function dodecahedron(s){
    const g=(1+Math.sqrt(5))/2;
    return [[s,s,s],[s,s,-s],[s,-s,s],[-s,s,s],[-s,-s,s],[-s,s,-s],[s,-s,-s],[-s,-s,-s],
            [0,s*g,s/g],[0,-s*g,s/g],[0,s*g,-s/g],[0,-s*g,-s/g],
            [s/g,0,s*g],[-s/g,0,s*g],[s/g,0,-s*g],[-s/g,0,-s*g],
            [s*g,s/g,0],[-s*g,s/g,0],[-s*g,-s/g,0],[s*g,-s/g,0]];
}
function segmentate(v,arr){
    const res=[];
    for(let k=0; k<arr.length; ++k){
        res.push(...geo.geodesic(v[arr[k][0]-1],v[arr[k][1]-1],12,true));
    }
    return res;
}
//0.1135
const H=dodecahedron(geo.radius*0.1134745);//spherical geometry geo.radius*0.1134745
//0.4275 in hyperbolic
const base={
    center:[0,0,0],
    segment:segmentate(H,[[1,17],[1,13],[1,9],[2,17],[2,11],[2,15],[3,10],[3,13],[3,20],[4,9],[4,14],[4,18],[5,10],[5,14],[5,19],[6,11],[6,18],[6,16],[7,15],[7,12],[7,20],[8,12],[8,16],[8,19],[9,11],[10,12],[19,18],[13,14],[15,16],[17,20]]),
    normal:[geo.midpoint([H[0],H[1],H[16],H[8],H[10]]),
            geo.midpoint([H[0],H[3],H[8],H[12],H[13]]),
            geo.midpoint([H[2],H[4],H[9],H[12],H[13]]),
            geo.midpoint([H[3],H[4],H[17],H[18],H[13]]),
            geo.midpoint([H[3],H[5],H[8],H[10],H[17]]),
            geo.midpoint([H[5],H[7],H[15],H[17],H[18]]),
            vectorneg(geo.midpoint([H[0],H[1],H[16],H[8],H[10]])),
            vectorneg(geo.midpoint([H[0],H[3],H[8],H[12],H[13]])),
            vectorneg(geo.midpoint([H[2],H[4],H[9],H[12],H[13]])),
            vectorneg(geo.midpoint([H[3],H[4],H[17],H[18],H[13]])),
            vectorneg(geo.midpoint([H[3],H[5],H[8],H[10],H[17]])),
            vectorneg(geo.midpoint([H[5],H[7],H[15],H[17],H[18]]))]
    };
function instantiate(V){
    const res=[];
    for(const v of V){
        res.push(v[0],v[1],v[2]);
    }
    return res;
}
const poly=[];
function tessellation(level){
    let stk=[base];
    let pdata=[base];
    function create(){
        const newdata=[];
        for(const p of pdata){
            //pdataの周りで。
            for(let k=0; k<p.normal.length; ++k){
                const n=p.normal[k];
                let safe=true;
                const cent=geo.refrection(p.center,n);
                if(stk.findIndex(e=>geo.distance(e.center,cent)<geo.radius*0.3)!=-1){
                    // safe=false;
                }
                if(safe){
                    newdata.push({
                        center:cent,
                        segment:geo.refrectionGroup(p.segment,n),
                        normal:geo.refrectionGroup(p.normal,n)
                    });
                }
            }
        }
        stk.push(...newdata);
        pdata=newdata;
    }
    for(let k=0; k<level-1; ++k){
    create();
    }
    for(const s of stk){
        poly.push(...s.segment);
    }
    // alert(stk.length);
}
tessellation(6);
vert=instantiate(poly);
wgpu.bindvertex(vert);

main();
