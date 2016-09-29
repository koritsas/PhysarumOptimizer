import edu.koritsas.slimemold.shapefile.IrrigationNetwork;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;

import java.io.IOException;
import java.util.List;

/**
 * Created by ilias on 29/9/2016.
 */
public class Main {
    public static void main(String[] args) {

        IrrigationNetwork network = new IrrigationNetwork("C:/Users/ilias/Desktop/SlimeTest/H.shp","C:/Users/ilias/Desktop/SlimeTest/WS.shp","C:/Users/ilias/Desktop/SlimeTest/P.shp");

        Graph graph=null;
        try {
            graph =network.getGraph();
        } catch (IOException e) {
            e.printStackTrace();
        }

        List<Node> sourceNodes =network.getWaterSource();
        List<Node> sinkNodes = network.getHydrants();
        sinkNodes=sinkNodes.subList(5,6);

        SlimeMold slimeMold = new SlimeMold(graph,sourceNodes,sinkNodes,2,1.8,500);
        slimeMold.execute();

        slimeMold.showFlowDiagram();
        slimeMold.showConductivityMap();

    }
}
