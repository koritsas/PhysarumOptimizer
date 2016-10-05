package edu.koritsas.slimemold;

import com.vividsolutions.jts.geom.Geometry;
import edu.koritsas.slimemold.shortestpath.PhysarumPolycephalumSP;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.graph.structure.Node;
import org.opengis.feature.simple.SimpleFeature;

/**
 * Created by ilias on 5/10/2016.
 */
public class SlimeMoldSP extends PhysarumPolycephalumSP {
    public SlimeMoldSP(Graph graph, Node sourceNode, Node sinkNode, double Io, double γ, int numberOfIterations) {
        super(graph, sourceNode, sinkNode, Io, γ, numberOfIterations);
    }

    @Override
    public double getEdgeCost(Edge e) {
        SimpleFeature f = (SimpleFeature) e.getObject();
        Geometry geom = (Geometry) f.getDefaultGeometry();

        return geom.getLength();
    }
}
