namespace Simulation.FEM {
    public class GlobalData {
        public double GlobalLengthGrid; //L szerokoœæ siatki
        public double GlobalHeightGrid; //H d³ugoœæ siatki 

        public int NHeight; //nH liczba wez³ów na wysokoœci
        public int NLength; //nL liczba wêz³ów na szerokoœci 

        public int NodesCount => NHeight * NLength;
        public int ElementsCount => (NHeight - 1) * (NLength - 1);
    }
}