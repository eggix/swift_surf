import ddf.minim.analysis.*;
//import ddf.minim.*;
import peasy.*;
import queasycam.*;

PeasyCam pcam;
QueasyCam qcam;
PFont mono;


boolean cam = false;
boolean game=true;

final int  N = 128;
float[][] field = new float[N][N];
float[][] fieldImg = new float[N][N];
float[][] fourierReal = new float[N][N];
float[][] fourierImg = new float[N][N];
float kcrit = 1;//N/100.*dk;
float dt = 0.01;
float eps = 1;
float alpha = 0;

float L = 10.0*(2.0*PI/kcrit);
float dx = L/N;
float dk = 2.0*PI/L;

final int ResX=1280;
final int ResY=720;
final float m=70;
final float E=1.5;
float max,min;
float score=0;
boolean gameOver=false;
boolean surfacelock=true;
//PShape goalPic=loadShape("GOAL.svg");

final float creep = 0.07;
final float penalty=20;
final int overSample=4;
final int colorMax=200;
final int colorMin=100;
final int backgroundColor=0;
final int saturation=200;
final int brightness=200;
final float range = 64;
final float power =40;
final float scale = 10;
final float vMax=1;
final int gameLength=5;
int worldCount=0;
FFT fft = new FFT(N,dk);

void fourierForward(){
  float[][] tempReal = new float[N][N];
  float[][] tempImg = new float[N][N];
  for(int i = 0; i < N; i++){
    fft.forward(field[i]);
    for(int k = 0; k < N; k++){
      //tempReal[k][i] =  pow(-1,k+i)*fft.getSpectrumReal()[(k+N/2)%N];
      //tempImg[k][i] =  pow(-1,k+i)*fft.getSpectrumImaginary()[(k+N/2)%N];
      tempReal[k][i] =  fft.getSpectrumReal()[k];
      tempImg[k][i] =  fft.getSpectrumImaginary()[k];
    }
  }
  for(int i = 0; i < N; i++){
    fft.forward(tempReal[i],tempImg[i]);
    for(int k = 0; k < N; k++){
      //fourierReal[i][k] = pow(-1,k+i)*fft.getSpectrumReal()[(k+N/2)%N];
      //fourierImg[i][k] = pow(-1,k+i)*fft.getSpectrumImaginary()[(k+N/2)%N];
      fourierReal[i][k] = fft.getSpectrumReal()[k]/N;
      fourierImg[i][k] = fft.getSpectrumImaginary()[k]/N;   
    } 
  }
}

void fourierInverse(){
 float[][] tempReal = new float[N][N];
  float[][] tempImg = new float[N][N];
  for(int i = 0; i < N; i++){
    fft.forward(fourierReal[i],fourierImg[i]);
    for(int k = 0; k < N; k++){
      //tempReal[k][i] = pow(-1,k+i)*fft.getSpectrumReal()[(k+N/2)%N];
      //tempImg[k][i] = pow(-1,k+i)*fft.getSpectrumImaginary()[(k+N/2)%N];
      tempReal[k][i] =  fft.getSpectrumReal()[k];
      tempImg[k][i] =  fft.getSpectrumImaginary()[k];
    }
  }
  for(int i = 0; i < N; i++){
    fft.forward(tempReal[i],tempImg[i]);
    for(int k = 0; k < N; k++){
      //field[i][k] = pow(-1,k+i)*fft.getSpectrumReal()[(k+N/2)%N];
      //fieldImg[i][k] = pow(-1,k+i)*fft.getSpectrumImaginary()[(k+N/2)%N];
      field[i][k] = fft.getSpectrumReal()[k]/N;
      fieldImg[i][k] = fft.getSpectrumImaginary()[k]/N;
    } 
  }
}

void linprop(){
  for(int i = 0; i < N; i++){
    for(int k = 0; k < N; k++){
      float laplace = kcrit*kcrit - (-N/2+i)*dk*(-N/2+i)*dk - (-N/2+k)*dk*(-N/2+k)*dk;
      fourierReal[(i+N/2)%N][(k+N/2)%N] *= 1/(1 - dt*(eps - laplace*laplace));
      fourierImg[(i+N/2)%N][(k+N/2)%N] *= 1/(1 - dt*(eps - laplace*laplace));
      //fourierReal[i][k] *= 1/(1 - dt*(eps - laplace*laplace));
      //fourierImg[i][k] *= 1/(1 - dt*(eps - laplace*laplace));

    } 
  }
}

void swift(){
  fourierForward();
  linprop();
  fourierInverse();
  for(int i = 0; i < N; i++){
    for(int k = 0; k < N; k++){
      float value = field[i][k];
      field[i][k] += dt * (-value*value*value + alpha*value*value);
    }
  }
}

void init(){
  for (int x = 0; x < N; x++){
    for(int y = 0; y < N; y++){
      //field[x][y] = exp(-(x-N/2)*(x-N/2)/30.-(y-N/2)*(y-N/2)/10.);
     //field[x][y] = sin(3*2*PI*x/float(N))+sin(7*2*PI*y/float(N));
      //field[x][y] = random(-0.1,0.1);
      field[x][y] = noise(x,y);
      //field[x][y] = x*y/(N);
      //field[x][y]=0;
     }
  }
}



void setup(){
  //size(1280,720,P3D);
  fullScreen(P3D);
  colorMode(HSB);
  mono = loadFont("URWGothic-Book-48.vlw");
  textFont(mono);
  textSize(12);
  textAlign(CENTER, CENTER);
  init();
  
  strokeWeight(5);
      
  makeLine();
  if(cam)
    pcam = new PeasyCam(this, 010);
  else{
    qcam=new QueasyCam(this);
    reset(true);
    //cam=new PeasyCam(this,100);
    //qcam.position=new PVector(N/2,-field[N/2][N-1],N-1);
    //qcam.controllable=true;
    //qcam.pan=-PI/2;
    //qcam.tilt=-PI/6;
  }
}

void draw(){

   background(backgroundColor);
   //text("d²/dt² psi = c² d²/dr² psi",0,100,-3*N);
//  rotateX(PI/2);
  //println(frameRate);
 // println(qcam.getForward().x + "    " + qcam.getForward().y +"     "+qcam.getForward().z);
// println(field[0][0]);
  for(int i=0;i<overSample;i++)
    swift();
   alpha=0.5+0.5*sin(frameCount/60./2.0);
   //eps=1+0.4*sin(frameCount/60./3.3);
  
  makeLine();
  if(!cam)  adjustQCam();
  if(game)  gameStuff();
 
}

void gameStuff(){
  
  showBeacon();
  if(gameOver){
    background(0);
    text("Game Over! \n Your Score is: " + score + "\n Press \"l\" to reinitialize the field \n Press \"r\" to restart ",N/2,-20,0);
  }else{
    score-=qcam.velocity.z*(field[int(qcam.position.x)][int(qcam.position.z)]-min);
    
    text("score: " + round(score),N/2,-50,-2*N+qcam.position.z+100);
    if(worldCount==gameLength){
      text("GOAL",N/2,-10,0);
      //shape(goalPic,0,-40,0);
    }
    
    addScore();
    if(frameCount%(30*3) == 0)
      beacons.add(new beacon());
    if(worldCount>gameLength)
    gameOver=true;
  }
}

void adjustQCam(){
 if(!game && (qcam.position.x< 0 || qcam.position.z < 0 || qcam.position.x>=N-0.5 || qcam.position.z>=N-0.5))
    reset(false);
 if(game && (qcam.position.x< 0 || qcam.position.z > N-0.5 || qcam.position.x>=N-0.5 )){
   gameOver=true;
   return;
 }
 if(game && qcam.position.z < 0){
   qcam.position.z=N-1;
   worldCount++;
 }
 if(surfacelock){
   if(abs(min-max)>1e-3)
     qcam.speed=sqrt(2*abs(field[int(qcam.position.x)][int(qcam.position.z)]-max+creep)/abs(min-max)*vMax/m);
   else
     qcam.speed=creep;
    qcam.position=new PVector(qcam.position.x,-field[round(qcam.position.x)][round(qcam.position.z)]-1,qcam.position.z);
 }
}

void makeLine(){
  max = -1000;
  min =  1000;
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      if(field[i][j]>max)
        max=field[i][j];
       if(field[i][j]<min)
         min=field[i][j];
    }
  }
  for(int i=0;i<N-1;i++){
    for(int j=0;j<N-1;j++){
      stroke(color((0.5*field[i][j]+0.5*field[i+1][j]-max)/(min-max)*(colorMax-colorMin),saturation,brightness));
      line(i,-field[i][j],j,i+1,-field[i+1][j],j);
      if(game)  line(i,-field[i][j],j-N,i+1,-field[i+1][j],j-N);
      stroke(color((0.5*field[i][j]+0.5*field[i][j+1]-max)/(min-max)*(colorMax-colorMin),saturation,brightness));
      line(i,-field[i][j],j,i,-field[i][j+1],j+1);
      if(game)  line(i,-field[i][j],j-N,i,-field[i][j+1],j+1-N);
    }
    
  }
  for(int i=0;i<N-1;i++){
    stroke(color((0.5*field[i][N-1]+0.5*field[i+1][N-1]-max)/(min-max)*(colorMax-colorMin),saturation,brightness));
    line(i,-field[i][N-1],N-1,i+1,-field[i+1][N-1],N-1);
    if(game){
      line(i,-field[i][N-1],-1,i+1,-field[i+1][N-1],-1);
      stroke(color((0.5*field[i][0]+0.5*field[i][N-1]-max)/(min-max)*(colorMax-colorMin),saturation,brightness));
      line(i,-field[i][0],0,i,-field[i][N-1],-1);
    }
    stroke(color((0.5*field[N-1][i]+0.5*field[N-1][i]-max)/(min-max)*(colorMax-colorMin),saturation,brightness));
    line(N-1,-field[N-1][i],i,N-1,-field[N-1][i+1],i+1);
  }
}

void reset(boolean setup){
  if(setup){
    qcam.pan=-PI/2;
    
  }
  qcam.position=new PVector(N/2,-field[N/2][N-1],N-1);
  qcam.tilt=-PI/12;
  for(int i=0; i<5; i++)
      beacons.add(new beacon());
  
}

void interaction(){
  float y = qcam.position.z+10*qcam.getForward().z;
  float x = qcam.position.x+10*qcam.getForward().x;
  for(int i=0;i<range;i++){
    
    for(int j=0;j<range;j++){
      
      field[(round(x+i-range/2.)+N)%N][(round(y+j-range/2.)+N)%N]-=scale*exp(-((i-range/2.)*(i-range/2.)+(j-range/2.)*(j-range/2.))/power);//*(sin(sqrt((x-i+range/2.)*(x-i+range/2.)+(y-j+range/2.)*(y-j+range/2.))*2*PI*kcrit));
    }
  }
}

void keyPressed(){
  if(key==' '){
    interaction();
    score-=penalty;
  }
   if(key=='f'){
     //cam=!cam;
     surfacelock=!surfacelock;
     qcam.speed=1;
   }   
  if(key=='r'&&gameOver){
    score=0;
    gameOver=false;
    worldCount=0;
    reset(true);
  }
  if(key=='l'&&gameOver){
    init();
  }
}
void mousePressed(){
   if(LEFT==mouseButton)
     alpha*=1.1;
   if(RIGHT==mouseButton)
     alpha/=1.1; 
}