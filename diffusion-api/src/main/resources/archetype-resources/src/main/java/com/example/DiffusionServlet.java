package com.example;

import javax.servlet.annotation.WebServlet;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import com.google.gson.Gson;
import java.util.HashMap;
import java.util.Map;

@WebServlet(name = "DiffusionServlet", value = "/calculate_diffusion")
public class DiffusionServlet extends HttpServlet {

    private final DiffusionCalculator calculator = new DiffusionCalculator();
    private final Gson gson = new Gson();

    @Override
    protected void doPost(HttpServletRequest request, HttpServletResponse response) throws IOException {
        String temperatureStr = request.getParameter("temperature");
        String pressureStr = request.getParameter("pressure");
        String viscosityStr = request.getParameter("viscosity");

        response.setContentType("application/json");
        PrintWriter out = response.getWriter();

        try {
            if (temperatureStr == null || pressureStr == null || viscosityStr == null) {
                response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
                Map<String, String> error = new HashMap<>();
                error.put("error", "Missing parameters: temperature, pressure, or viscosity are required.");
                out.print(gson.toJson(error));
                return;
            }

            double temperature = Double.parseDouble(temperatureStr);
            double pressure = Double.parseDouble(pressureStr);
            double viscosity = Double.parseDouble(viscosityStr);

            double diffusionCoefficient = calculator.calculateDiffusionCoefficient(temperature, pressure, viscosity);

            Map<String, Double> result = new HashMap<>();
            result.put("diffusionCoefficient", diffusionCoefficient);
            out.print(gson.toJson(result));

        } catch (NumberFormatException e) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            Map<String, String> error = new HashMap<>();
            error.put("error", "Invalid parameter format: temperature, pressure, and viscosity must be numbers.");
            out.print(gson.toJson(error));
        } catch (IllegalArgumentException e) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            Map<String, String> error = new HashMap<>();
            error.put("error", e.getMessage());
            out.print(gson.toJson(error));
        } finally {
            out.close();
        }
    }
}