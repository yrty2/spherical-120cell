const geo=new geometry(1,3);
const wgpu=new WGPU(geo,8,`
struct Uniforms {
    constants:vec4<f32>,
    rot:vec4<f32>
}
@binding(0) @group(0) var<uniform> uni : Uniforms;

struct VertexOutput {
  @builtin(position) Position:vec4<f32>,
  @location(0) fragColor:vec3<f32>
}
fn qconj(p:vec4<f32>)->vec4<f32>{
  return vec4<f32>(p.x,-p.yzw);
}
fn qmul(p:vec4<f32>,q:vec4<f32>)->vec4<f32>{
  return p.x*q+vec4<f32>(-dot(p.yzw,q.yzw),cross(p.yzw,q.yzw)+q.x*p.yzw);
}
fn rotate(p:vec3<f32>,rot:vec4<f32>)->vec3<f32>{
  return qmul(qmul(rot,vec4<f32>(0,p)),qconj(rot)).yzw;
}
fn sphericaldist(p:vec3<f32>)->f32{
  return 2*atan(length(p));
}
fn klein(p:vec3<f32>)->vec3<f32>{
  return 2*p/(1+dot(p,p));
}
fn up(p:vec3<f32>)->vec3<f32>{
  return vec3<f32>(p.xy,sqrt(1-dot(p,p)))/(1-p.z)-vec3<f32>(0,0,1);
}
@vertex
fn main(@location(0) position:vec3<f32>)->VertexOutput{
  var output : VertexOutput;
  var p=rotate((position)+uni.constants.yzw,uni.rot);
  let dist:f32=sphericaldist(position);
  p=vec3<f32>(p.xy/(1.4142*p.z),p.z*0.00001);
  p.y*=uni.constants.x;
  if(p.z>0){
  output.Position=vec4<f32>(p,1);
  output.fragColor=vec3<f32>(1,1,1)/dist;
  }
  return output;
}
@fragment
fn fragmain(@location(0) fragColor:vec3<f32>)->@location(0) vec4<f32>{
  return vec4<f32>(fragColor,1);
}
`,"segment");
async function main(){
    await wgpu.initialize(canvas,["float32x3"]);
    function gameloop(){
        wgpu.background=[0,0,0];
        wgpu.uniform=[screen.width/screen.height,...cam,...rot];
        wgpu.render();
        requestAnimationFrame(gameloop);
    }
    gameloop();

}
