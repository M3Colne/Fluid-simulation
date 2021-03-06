/****************************************************************************************** 
 *	Chili DirectX Framework Version 16.07.20											  *	
 *	Game.cpp																			  *
 *	Copyright 2016 PlanetChili.net <http://www.planetchili.net>							  *
 *																						  *
 *	This file is part of The Chili DirectX Framework.									  *
 *																						  *
 *	The Chili DirectX Framework is free software: you can redistribute it and/or modify	  *
 *	it under the terms of the GNU General Public License as published by				  *
 *	the Free Software Foundation, either version 3 of the License, or					  *
 *	(at your option) any later version.													  *
 *																						  *
 *	The Chili DirectX Framework is distributed in the hope that it will be useful,		  *
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of						  *
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the						  *
 *	GNU General Public License for more details.										  *
 *																						  *
 *	You should have received a copy of the GNU General Public License					  *
 *	along with The Chili DirectX Framework.  If not, see <http://www.gnu.org/licenses/>.  *
 ******************************************************************************************/
#include "MainWindow.h"
#include "Game.h"
#include <cassert>

Game::Game( MainWindow& wnd )
	:
	wnd( wnd ),
	gfx( wnd )
{
	//Initialization
	prev_mousePos = wnd.mouse.GetPos();
	const int totalCells = N * N;
	density = new float[totalCells];
	prev_density = new float[totalCells];
	velocity = new Vec2[totalCells];
	prev_velocity = new Vec2[totalCells];

	for (int i = 0; i < totalCells; i++)
	{
		density[i] = 0.0f;
		prev_density[i] = 0.0f;
		velocity[i] = Vec2(0.0f, 0.0f);
		prev_velocity[i] = Vec2(0.0f, 0.0f);
	}

	//Using an initial scalar field to test
	/*for (int j = 1; j < N-1; j++)
	{
		for (int i = 1; i < N-1; i++)
		{
			const float x = float(i - int(N / 2));
			const float y = float(j - int(N / 2));

			density[GetId(i, j)] = pow(2.0f, -1.0f * Vec2(x, y).GetLengthSq());
		}
	}*/

	//Using an initial vector field to test
	//for (int j = 0; j < N; j++)
	//{
	//	for (int i = 0; i < N; i++)
	//	{
	//		const float x = float(i - int(N / 2));
	//		const float y = float(j - int(N / 2));

	//		//Arbitrary 2D function
	//		velocity[GetId(i, j)] = Vec2(x, y);
	//	}
	//}
}

Game::~Game()
{
	delete[] density;
	density = nullptr;
	delete[] prev_density;
	prev_density = nullptr;
	delete[] velocity;
	velocity = nullptr;
	delete[] prev_velocity;
	prev_velocity = nullptr;
}

void Game::Go()
{
	gfx.BeginFrame();	
	UpdateModel();
	ComposeFrame();
	gfx.EndFrame();
}

void Game::UpdateModel()
{
	//Delta time
	const float DT = ft.Mark();
	//Delta time

	//Reset
	if (wnd.kbd.KeyIsPressed('R'))
	{
		for (int i = 0; i < N * N; i++)
		{
			density[i] = 0.0f;
			velocity[i].x = 0.0f;
			velocity[i].y = 0.0f;
		}
	}

	//Pause
	if (wnd.kbd.KeyIsPressed('P'))
	{
		if (pauseInhib)
		{
			pause = !pause;
			pauseInhib = false;
		}
	}
	else
	{
		pauseInhib = true;
	}

	//Navier-Stoke equations for mass and velocity
	if (!pause)
	{
		//This is the material derivative
		DensitySolver(DPS, brushRadius, diffusionRate, DT);
		//This is the Navier-Stokes equation for velocity
		VelocitySolver(velocityScalar, brushRadius, viscosityRate, DT);

		//This is to prevent having too small values like 0.0001 for density and velocity
		for (int i = 0; i < N * N; i++)
		{
			const float d = density[i];
			density[i] = d < minDensity ? 0.0f : d;
			velocity[i] = velocity[i].GetLengthSq() < minVelocitySq ? Vec2(0.0f, 0.0f) : velocity[i];
		}
	}

	//Updating the mouse position
	prev_mousePos = wnd.mouse.GetPos();
}

int Game::GetId(int i, int j)
{
	return j * N + i;
}

float Game::LinearInterpolation(float a, float b, float x)
{
	return (1.0f-x)*a + x*b;
}

void Game::DrawDensity()
{
	//Drawing every cell
	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			const int A = int(j * cellDimension);
			const int B = int(i * cellDimension);
			for (int J = A; J < A+cellDimension; J++)
			{
				for (int I = B; I < B+cellDimension; I++)
				{
					const int id = GetId(i, j);
					int greyScale = int(255 * density[id]);
					if (greyScale > 255)
					{
						greyScale = 255;
					}
					gfx.PutPixel(I, J, greyScale, greyScale, greyScale);
				}
			}
		}
	}
}

void Game::DrawVelocities()
{
	//The velocities colors are relative to the maximum length
	const float maxLengthDrawn = cellDimension / 2.0f;
	const float maxLengthDrawnSq = maxLengthDrawn * maxLengthDrawn;

	for (int j = 0; j < N; j++)
	{
		for (int i = 0; i < N; i++)
		{
			const Vec2 center(cellDimension * i + maxLengthDrawn, cellDimension * j + maxLengthDrawn);

			const float vel = velocity[GetId(i, j)].GetLengthSq();
			Color color;
			if (vel > maxLengthDrawn * maxLengthDrawn)
			{
				color = Colors::Magenta;
			}
			else if (vel > maxLengthDrawnSq * 16.0f / 25.0f)
			{
				color = Colors::Red;
			}
			else if (vel > maxLengthDrawnSq * 9.0f / 25.0f)
			{
				color = Colors::Yellow;
			}
			else if (vel > maxLengthDrawnSq * 4.0f / 25.0f)
			{
				color = Colors::Green;
			}
			else if (vel > maxLengthDrawnSq / 25.0f)
			{
				color = Colors::Cyan;
			}
			else
			{
				color = Colors::Blue;
			}
			gfx.DrawLine(center, center + velocity[GetId(i, j)].GetNormalizedTo(maxLengthDrawn), color);
		}
	}
}

void Game::DensitySolver(float brushAmountPerSec, float brushRadius, float diffRate, float dt)
{
	AddDensity(brushAmountPerSec, brushRadius + 0.5f, dt);
	std::swap(density, prev_density);
	Diffusion(diffRate, dt);
	std::swap(density, prev_density);
	Advection(dt);
}

void Game::VelocitySolver(float scalar, float brushRadius, float viscRate, float dt)
{
	AddVelocity(scalar, brushRadius + 0.5f);
	std::swap(velocity, prev_velocity);
	Viscosity(viscRate, dt);
	////Maybe another project comes here because it says it may work better
	std::swap(velocity, prev_velocity);
	Convection(dt);
	PressureProjection();
}

void Game::DensityBoundaryCondition()
{
	for (int j = 1; j < N - 1; j++)
	{
		density[GetId(0, j)] = density[GetId(1, j)];
		density[GetId(N - 1, j)] = density[GetId(N - 2, j)];
	}
	for (int i = 1; i < N - 1; i++)
	{
		density[GetId(i, 0)] = density[GetId(i, 1)];
		density[GetId(i, N - 1)] = density[GetId(i, N - 2)];
	}
	density[GetId(0, 0)] = (density[GetId(1, 0)] + density[GetId(0, 1)]) / 2.0f;
	density[GetId(N-1, 0)] = (density[GetId(N-2, 0)] + density[GetId(N-1, 1)]) / 2.0f;
	density[GetId(0, N-1)] = (density[GetId(0, N-2)] + density[GetId(1, N-1)]) / 2.0f;
	density[GetId(N-1, N-1)] = (density[GetId(N-2, N-1)] + density[GetId(N-1, N-2)]) / 2.0f;
}

void Game::VelocityBoundaryCondition()
{
	for (int j = 1; j < N - 1; j++)
	{
		velocity[GetId(0, j)].x = velocity[GetId(1, j)].x * -1;
		velocity[GetId(N - 1, j)].x = velocity[GetId(N - 2, j)].x * -1;
	}
	for (int i = 1; i < N - 1; i++)
	{
		velocity[GetId(i, 0)].y = velocity[GetId(i, 1)].y * -1;
		velocity[GetId(i, N - 1)].y = velocity[GetId(i, N - 2)].y * -1;
	}
	velocity[GetId(0, 0)] = (velocity[GetId(1, 0)] + velocity[GetId(0, 1)]) / 2.0f;
	velocity[GetId(N - 1, 0)] = (velocity[GetId(N - 2, 0)] + velocity[GetId(N - 1, 1)]) / 2.0f;
	velocity[GetId(0, N - 1)] = (velocity[GetId(0, N - 2)] + velocity[GetId(1, N - 1)]) / 2.0f;
	velocity[GetId(N - 1, N - 1)] = (velocity[GetId(N - 2, N - 1)] + velocity[GetId(N - 1, N - 2)]) / 2.0f;
}

void Game::AddDensity(float AmountPerSec, float radius, float dt)
{
	if (wnd.mouse.LeftIsPressed())
	{
		Vec2 mousePos(float(wnd.mouse.GetPosX()), float(wnd.mouse.GetPosY()));
		mousePos /= cellDimension;

		const float radiusSq = radius * radius;
		int xStart = int(mousePos.x - radius);
		int yStart = int(mousePos.y - radius);
		int xEnd = int(mousePos.x + radius);
		int yEnd = int(mousePos.y + radius);
		if (xStart < 1)
		{
			xStart = 1;
		}
		if (yStart < 1)
		{
			yStart = 1;
		}
		if (xEnd > N - 1)
		{
			xEnd = N - 1;
		}
		if (yEnd > N - 1)
		{
			yEnd = N - 1;
		}
		
		for (int j = yStart; j < yEnd; j++)
		{
			for (int i = xStart; i < xEnd; i++)
			{
				if ((mousePos.x - i) * (mousePos.x - i) + (mousePos.y - j) * (mousePos.y - j) <= radiusSq)
				{
					density[GetId(i, j)] += dt * AmountPerSec;
				}
			}
		}
	}
}

void Game::AddVelocity(float scalar, float radius)
{
	if (wnd.mouse.RightIsPressed())
	{
		Vec2 mousePos(float(wnd.mouse.GetPosX()), float(wnd.mouse.GetPosY()));
		mousePos /= cellDimension;

		const float radiusSq = radius * radius;
		int xStart = int(mousePos.x - radius);
		int yStart = int(mousePos.y - radius);
		int xEnd = int(mousePos.x + radius);
		int yEnd = int(mousePos.y + radius);
		if (xStart < 1)
		{
			xStart = 1;
		}
		if (yStart < 1)
		{
			yStart = 1;
		}
		if (xEnd > N - 1)
		{
			xEnd = N - 1;
		}
		if (yEnd > N - 1)
		{
			yEnd = N - 1;
		}

		Vec2 mouseTimeDiff = mousePos - Vec2(float(prev_mousePos.x), float(prev_mousePos.y))/cellDimension;

		for (int j = yStart; j < yEnd; j++)
		{
			for (int i = xStart; i < xEnd; i++)
			{
				if ((mousePos.x - i) * (mousePos.x - i) + (mousePos.y - j) * (mousePos.y - j) <= radiusSq)
				{
					velocity[GetId(i, j)] += mouseTimeDiff * scalar;
				}
			}
		}
	}
}

void Game::Diffusion(float diffRate, float dt)
{
	const float a = dt * diffRate;
	for (int iteration = 0; iteration < 20; iteration++)
	{
		for (int j = 1; j <= n; j++)
		{
			for (int i = 1; i <= n; i++)
			{
				const int id = GetId(i, j);
				density[id] = (prev_density[id] + a *
					(density[id - 1] + density[id + 1] + density[id + N] + density[id - N]))/(1+4*a);
			}
		}
		DensityBoundaryCondition();
	}
}

void Game::Viscosity(float viscosityRate, float dt)
{
	const float a = dt * viscosityRate;
	for (int iteration = 0; iteration < 20; iteration++)
	{
		for (int j = 1; j <= n; j++)
		{
			for (int i = 1; i <= n; i++)
			{
				const int id = GetId(i, j);
				velocity[id] = (prev_velocity[id] +	
					(velocity[id - 1] + velocity[id + 1] + velocity[id + N] + velocity[id - N]) * a) / (1 + 4 * a);
			}
		}
		VelocityBoundaryCondition();
	}
}

void Game::Advection(float dt)
{
	//Simple Backwards Euler
	//Semi-Lagrangian advection (going backwards in time)
	for (float j = 1.0f; j <= n; j++)
	{
		for (float i = 1.0f; i <= n; i++)
		{
			//Going backward in time
			Vec2 pos(i, j); //Position in a 1-n * 1-n grid
							//Also, it should be pos(float(i+0.5f), float(j+0.5f)) but for some reason it's buggy
			pos -= velocity[GetId(int(i), int(j))] * dt;

			//Constraints
			if (pos.x < 0.0f)
			{
				pos.x = 0.0f;
			}
			else if (pos.x > N)
			{
				pos.x = float(N);
			}
			if (pos.y < 0.0f)
			{
				pos.y = 0.0f;
			}
			else if (pos.y > N)
			{
				pos.y = float(N);
			}

			//Interpolating the particle density around his 4 neighbors
			const int nPosX = int(pos.x);
			const int nPosY = int(pos.y);
			const float fracY = pos.y - nPosY;
			const float Y1 = LinearInterpolation(prev_density[GetId(nPosX, nPosY)], prev_density[GetId(nPosX, nPosY + 1)], fracY);
			const float Y2 = LinearInterpolation(prev_density[GetId(nPosX+1, nPosY)], prev_density[GetId(nPosX+1, nPosY + 1)], fracY);
			density[GetId(int(i), int(j))] = LinearInterpolation(Y1, Y2, pos.x - nPosX);
		}
	}
	DensityBoundaryCondition();
}

void Game::Convection(float dt)
{
	//Simple Backwards Euler
	//Semi-Lagrangian advection (going backwards in time)
	for (float j = 1; j <= n; j++)
	{
		for (float i = 1; i <= n; i++)
		{
			//Going backward in time
			const int id = GetId(int(i), int(j));
			Vec2 pos(i, j); //Position in a 0-N * 0-N grid
			pos -= velocity[id] * dt;

			//Constraints
			if (pos.x < 0.0f)
			{
				pos.x = 0.0f;
			}
			else if (pos.x > N)
			{
				pos.x = N;
			}
			if (pos.y < 0.0f)
			{
				pos.y = 0.0f;
			}
			else if (pos.y > N)
			{
				pos.y = N;
			}

			//Interpolating the particle density around his 4 neighbors
			const int nPosX = int(pos.x);
			const int nPosY = int(pos.y);
			const float fracY = pos.y - nPosY;
			const float Y1 = LinearInterpolation(prev_velocity[GetId(nPosX, nPosY)].x, prev_velocity[GetId(nPosX, nPosY + 1)].x, fracY);
			const float Y2 = LinearInterpolation(prev_velocity[GetId(nPosX + 1, nPosY)].x, prev_velocity[GetId(nPosX + 1, nPosY + 1)].x, fracY);
			velocity[id].x = LinearInterpolation(Y1, Y2, pos.x - nPosX);
			const float Y3 = LinearInterpolation(prev_velocity[GetId(nPosX, nPosY)].y, prev_velocity[GetId(nPosX, nPosY + 1)].y, fracY);
			const float Y4 = LinearInterpolation(prev_velocity[GetId(nPosX + 1, nPosY)].y, prev_velocity[GetId(nPosX + 1, nPosY + 1)].y, fracY);
			velocity[id].y = LinearInterpolation(Y3, Y4, pos.x - nPosX);
		}
	}
	VelocityBoundaryCondition();
}

void Game::PressureProjection()
{
	//Using Stam's pressure solver
	//I believe he is using an iterative solver called Gauss-Seidel relaxation to obtain the pressure

	//Observation:
	//The prev_velocity.y is supposed to hold the divergence values for each cell
	//The prev_velocity.x is supposed to hold the pressure values for each cell

	//Divergence and initializing the pressure
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= n; i++)
		{
			prev_velocity[GetId(i, j)].y = -0.5f * (velocity[GetId(i + 1, j)].x - velocity[GetId(i - 1, j)].x +
				velocity[GetId(i, j + 1)].y - velocity[GetId(i, j - 1)].y);
			prev_velocity[GetId(i, j)].x = 0.0f;
		}
	}

	//Divergence and pressure boundary condition
	{
		//Vertical boundaries
		for (int j = 0; j < N; j++)
		{
			//Divergence
			prev_velocity[GetId(0, j)].y = prev_velocity[GetId(1, j)].y;
			prev_velocity[GetId(N - 1, j)].y = prev_velocity[GetId(N - 2, j)].y;
			//Pressure
			prev_velocity[GetId(0, j)].x = 0.0f;
			prev_velocity[GetId(N - 1, j)].x = 0.0f;
		}
		//Horizontal boundaries
		for (int i = 0; i < N; i++)
		{
			//Divergence
			prev_velocity[GetId(i, 0)].y = prev_velocity[GetId(i, 1)].y;
			prev_velocity[GetId(i, N - 1)].y = prev_velocity[GetId(i, N - 2)].y;
			//Pressure
			prev_velocity[GetId(i, 0)].x = 0.0f;
			prev_velocity[GetId(i, N - 1)].x = 0.0f;
		}
		//Corners
		prev_velocity[GetId(0, 0)] = (prev_velocity[GetId(1, 0)] + prev_velocity[GetId(0, 1)]) / 2.0f;
		prev_velocity[GetId(N - 1, 0)] = (prev_velocity[GetId(N - 2, 0)] + prev_velocity[GetId(N - 1, 1)]) / 2.0f;
		prev_velocity[GetId(0, N - 1)] = (prev_velocity[GetId(0, N - 2)] + prev_velocity[GetId(1, N - 1)]) / 2.0f;
		prev_velocity[GetId(N - 1, N - 1)] = (prev_velocity[GetId(N - 2, N - 1)] + prev_velocity[GetId(N - 1, N - 2)]) / 2.0f;
	}

	//Pressure solving
	{
		for (int k = 0; k < 20; k++)
		{
			for (int j = 1; j <= n; j++)
			{
				for (int i = 1; i <= n; i++)
				{
					prev_velocity[GetId(i, j)].x = (prev_velocity[GetId(i, j)].y + prev_velocity[GetId(i + 1, j)].x + prev_velocity[GetId(i - 1, j)].x +
						prev_velocity[GetId(i, j + 1)].x + prev_velocity[GetId(i, j - 1)].x) / 4;
				}
			}

			//Pressure boundary condition
			//Vertical boundaries
			for (int j = 0; j < N; j++)
			{
				prev_velocity[GetId(0, j)].x = 0.0f;
				prev_velocity[GetId(N - 1, j)].x = 0.0f;
			}
			//Horizontal boundaries
			for (int i = 0; i < N; i++)
			{
				prev_velocity[GetId(i, 0)].x = 0.0f;
				prev_velocity[GetId(i, N - 1)].x = 0.0f;
			}
			//Corners
			prev_velocity[GetId(0, 0)].x = (prev_velocity[GetId(1, 0)].x + prev_velocity[GetId(0, 1)].x) / 2.0f;
			prev_velocity[GetId(N - 1, 0)].x = (prev_velocity[GetId(N - 2, 0)].x + prev_velocity[GetId(N - 1, 1)].x) / 2.0f;
			prev_velocity[GetId(0, N - 1)].x = (prev_velocity[GetId(0, N - 2)].x + prev_velocity[GetId(1, N - 1)].x) / 2.0f;
			prev_velocity[GetId(N - 1, N - 1)].x = (prev_velocity[GetId(N - 2, N - 1)].x + prev_velocity[GetId(N - 1, N - 2)].x) / 2.0f;
		}
	}

	//Using pressure to correct the velocities
	for (int j = 1; j <= n; j++)
	{
		for (int i = 1; i <= n; i++)
		{
			velocity[GetId(i, j)].x -= 0.5f * (prev_velocity[GetId(i + 1, j)].x - prev_velocity[GetId(i - 1, j)].x);
			velocity[GetId(i, j)].y -= 0.5f * (prev_velocity[GetId(i, j + 1)].x - prev_velocity[GetId(i, j - 1)].x);
		}
	}
	VelocityBoundaryCondition();
}

void Game::ComposeFrame()
{
	if (wnd.kbd.KeyIsPressed('D'))
	{
		if (drawDInhib)
		{
			drawD = !drawD;
			drawDInhib = false;
		}
	}
	else
	{
		drawDInhib = true;
	}
	if (drawD)
	{
		DrawDensity();
	}

	if (wnd.kbd.KeyIsPressed('V'))
	{
		if (drawVInhib)
		{
			drawV = !drawV;
			drawVInhib = false;
		}
	}
	else
	{
		drawVInhib = true;
	}
	if (drawV)
	{
		DrawVelocities();
	}
}
