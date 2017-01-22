ArrayList<beacon> beacons = new ArrayList();

public class beacon{
  PVector pos;
  color c;
  public beacon(){
       int x = int(random(0,N));
       int z = int(random(0,N));
       c = color(random(colorMin,colorMax),saturation,brightness); 
       pos = new PVector(x,-field[x][z],z);
  }
}
  public  void addScore(){
    beacon toRemove=null;
    for(beacon b: beacons){
      if(PVector.dist(qcam.position,b.pos) < 5){
         score += 20;
         toRemove=b;
      }
    }
    beacons.remove(toRemove);
  }
  
  public  void showBeacon(){
    for(beacon b: beacons){
      strokeWeight(8);
      stroke(b.c);
      line(b.pos.x,b.pos.y,b.pos.z,b.pos.x,b.pos.y-100,b.pos.z);
      line(b.pos.x,b.pos.y,b.pos.z-N,b.pos.x,b.pos.y-100,b.pos.z-N);

    }
  }