namespace Simulation.FEM {
    public class GlobalData {
        public double GlobalLengthGrid; //L szeroko�� siatki
        public double GlobalHeightGrid; //H d�ugo�� siatki 

        public int NHeight; //nH liczba wez��w na wysoko�ci
        public int NLength; //nL liczba w�z��w na szeroko�ci 

        public int NodesCount => NHeight * NLength;
        public int ElementsCount => (NHeight - 1) * (NLength - 1);
    }
}