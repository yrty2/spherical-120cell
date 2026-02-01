const keys={
    w:false,
    a:false,
    s:false,
    d:false
};
function keyframe(){
    let mv=[0,0,0];
    if(keys.w){
        mv[2]++;
    }
    if(keys.a){
        mv[0]--;
    }
    if(keys.s){
        mv[2]--;
    }
    if(keys.d){
        mv[0]++;
    }
    if(keys.shiftl){
        mv[1]--;
    }
    if(keys.space){
        mv[1]++;
    }
    const s=vectorlength(mv);
    if(s>0){
        vert=[];
        const t=q256.rotate(vectormul(mv,-1/(60*s)),q256.conj(rot));
        for(let k=0; k<poly.length; ++k){
            //if(geo.curvature>0){
            //poly[k]=geo.cliffordtranslate(poly[k],t);
            //}else{
                poly[k]=geo.translate(poly[k],t);
            //}
            wgpu.vertex[3*k]=poly[k][0];
            wgpu.vertex[3*k+1]=poly[k][1];
            wgpu.vertex[3*k+2]=poly[k][2];
        }
    }
    requestAnimationFrame(keyframe);
}
keyframe();
window.addEventListener("keydown",e=>{
    if(e.code=="KeyW"){
        keys.w=true;
    }
    if(e.code=="KeyA"){
        keys.a=true;
    }
    if(e.code=="KeyS"){
        keys.s=true;
    }
    if(e.code=="KeyD"){
        keys.d=true;
    }
    if(e.code=="Space"){
        keys.space=true;
    }
    if(e.code=="ShiftLeft"){
        keys.shiftl=true;
    }
});
window.addEventListener("keyup",e=>{
    if(e.code=="KeyW"){
        keys.w=false;
    }
    if(e.code=="KeyA"){
        keys.a=false;
    }
    if(e.code=="KeyS"){
        keys.s=false;
    }
    if(e.code=="KeyD"){
        keys.d=false;
    }
    if(e.code=="Space"){
        keys.space=false;
    }
    if(e.code=="ShiftLeft"){
        keys.shiftl=false;
    }
});
canvas.addEventListener("pointermove",e=>{
    canvas.requestPointerLock =
    canvas.requestPointerLock || canvas.mozRequestPointerLock;
    canvas.requestPointerLock();
    const dv=vectormul([e.movementX,e.movementY],1/1000);
    rot=q256.mul(q256.camrot(dv),rot);
});