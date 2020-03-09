namespace Simulation.FEM {
    public class Node {
        public double GlobalX, GlobalY;
        public double Value; //t
        public int X;
        public int Y;
        public int Id;
        public int BoundaryCondition;
        public bool IsBoundary => BoundaryCondition != 0;
        
    }
}