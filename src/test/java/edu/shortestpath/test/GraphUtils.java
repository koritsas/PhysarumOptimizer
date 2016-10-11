package edu.shortestpath.test;

import org.geotools.data.DefaultTransaction;
import org.geotools.data.FileDataStore;
import org.geotools.data.FileDataStoreFinder;
import org.geotools.data.Transaction;
import org.geotools.data.collection.ListFeatureCollection;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.graph.structure.Edge;
import org.geotools.graph.structure.Graph;
import org.geotools.map.FeatureLayer;
import org.geotools.map.Layer;
import org.geotools.map.MapContent;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.geotools.styling.SLD;
import org.geotools.styling.Style;
import org.geotools.swing.JMapFrame;
import org.opengis.feature.simple.SimpleFeature;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by ilias on 29/9/2016.
 */
public class GraphUtils {

    public static File createShpFromGraph(Graph graph,String path) throws IOException {

        List<Edge> edges = new ArrayList<Edge>(graph.getEdges());


        List<SimpleFeature> features = new ArrayList<SimpleFeature>();
        for (Edge e : edges) {
            SimpleFeature f = (SimpleFeature) e.getObject();
            features.add(f);
        }

        SimpleFeatureCollection featureCollection = new ListFeatureCollection(features.get(0).getFeatureType(), features);

        ShapefileDataStoreFactory factory = new ShapefileDataStoreFactory();

        File newFile = new File(path);

        ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();

        Map<String, Serializable> params = new HashMap<String, Serializable>();
        params.put("url", newFile.toURI().toURL());
        params.put("create spatial index", Boolean.TRUE);

        ShapefileDataStore newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);


        newDataStore.createSchema(features.get(0).getFeatureType());

        /*
         * You can comment out this line if you are using the createFeatureType method (at end of
         * class file) rather than DataUtilities.createType
         */
        newDataStore.forceSchemaCRS(DefaultGeographicCRS.WGS84);

        Transaction transaction = new DefaultTransaction("create");

        String typeName = newDataStore.getTypeNames()[0];
        SimpleFeatureSource featureSource = newDataStore.getFeatureSource(typeName);

        if (featureSource instanceof SimpleFeatureStore) {
            SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;

            featureStore.setTransaction(transaction);
            try {
                featureStore.addFeatures(featureCollection);
                transaction.commit();

            } catch (Exception problem) {
                problem.printStackTrace();
                transaction.rollback();

            } finally {
                transaction.close();
            }

        } else {
            System.out.println(typeName + " does not support read/write access");
            System.exit(1);
        }

        return newFile;
    }

    public static void readShapeFile(File file) throws IOException {




        FileDataStore store = FileDataStoreFinder.getDataStore(file);
        SimpleFeatureSource source = store.getFeatureSource();

        // Create a map content and add our shapefile to it
        MapContent map = new MapContent();
        map.setTitle("Quickstart");

        Style style = SLD.createSimpleStyle(source.getSchema());
        Layer layer = new FeatureLayer(source, style);
        map.addLayer(layer);

        // Now display the map
        JMapFrame.showMap(map);

    }

    public static void visualizeGraph(Graph graph) throws IOException {

        File file =File.createTempFile("Graph",".shp");

        file =createShpFromGraph(graph,file.getAbsolutePath());

        readShapeFile(file);

        file.deleteOnExit();

    }
}
